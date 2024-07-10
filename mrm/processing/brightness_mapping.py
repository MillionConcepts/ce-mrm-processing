from pathlib import Path
from typing import Literal, Sequence
import warnings

from astropy.io import fits
import cv2
from dustgoggles.func import gmap
import numpy as np
import pdr
import pyarrow as pa
from mrm.processing.hdus import to_scaled_hdu
from pyarrow import compute as pc
from pyarrow import parquet as pq
from scipy.interpolate import interpn

from mrm.shared.constants import CHANNEL_FREQ
from mrm.shared.fill_invalid import fill_invalid
from mrm.shared.pqutils import pq_metadict

# latitude cutoffs of input maps (they should all have full global longitude
# coverage, although not all the same center)
SOURCE_LAT = {
    'hparam': (70, -70),
    'slope': (90, -90),
    'azimuth': (-90, 90),
    'albedo': (70, -70)
}
# resolution in deg/pix of the input maps
SOURCE_DPP = {
    'hparam': 1 / 64, 'slope': 1 / 64, 'azimuth': 1 / 64, 'albedo': 1 / 64
}
# filenames of input maps
SRCMAP_FNS = {
    'hparam': 'H-parameter-64ppd.fits',
    'slope': 'SLDEM2015_64_SL_90S_90N_000_360.LBL',
    'azimuth': 'SLDEM2015_64_AZ_90S_90N_000_360.LBL',
    'albedo': 'wac_hapke_604nm_70N70S_64ppd.fits'
}
SRCTABLE_FN = 'tb_table.parquet'
# paths to source / intermediate maps, populated at runtime by make_filepaths()
SRC_PATHS, INT_PATHS = {}, {}
# type alias for map names
SourceMap = Literal["hparam", "slope", "azimuth", "albedo"]


def _scalefactor(mapname: SourceMap, dpp: float) -> int:
    if (scale := dpp / SOURCE_DPP[mapname]) % 1 != 0:
        raise ValueError(
            f"Requested resolution of {dpp} dpp does not scale evenly to the "
            f"{SOURCE_DPP[mapname]} resolution of the {mapname} map."
        )
    return int(scale)


# WARNING: expects that maparr has already been resized.
def _intermediate_hdul(
    maparr: np.ndarray, mapname: SourceMap, dpp: float
) -> fits.HDUList:
    lbounds = SOURCE_LAT[mapname]
    lat = np.arange(lbounds[0] - dpp, lbounds[1] - dpp, -dpp, dtype='f4')
    lon = np.arange(-180, 180, dpp, dtype='f4')
    # noinspection PyTypeChecker
    hdus = [
        fits.PrimaryHDU(),
        fits.ImageHDU(maparr, name=mapname),
        fits.ImageHDU(lat, name="LATITUDE"),
        fits.ImageHDU(lon, name="LONGITUDE")
    ]
    hdul = fits.HDUList(hdus)
    hdul[mapname].header['RESAMPLE'] = _scalefactor(mapname, dpp)
    hdul[mapname].header['dpp'] = int(1 / dpp)
    return hdul


def _shrink(
    big: np.ndarray, mapname: SourceMap, dpp: float
) -> np.ndarray:
    if (scl := _scalefactor(mapname, dpp)) == 1:
        return big
    small = np.zeros((big.shape[0] // scl, big.shape[1] // scl), np.float32)
    return cv2.resize(
        big, np.flip(small.shape), small, interpolation=cv2.INTER_AREA
    )


def _manipulate_hparam(sourcefile, dpp, valid_max, model_bounds):
    hparam = fits.open(sourcefile)[0].data.astype('f4')
    reflat = np.flip(np.arange(-90, 90, SOURCE_DPP['hparam']))
    # NOTE: AFAIK, +/- 70 lat is the valid range of the map. Someone
    #  just inset it in an array sized for -90 to 90. It would be nice to
    #  have metadata!
    lmax, lmin = SOURCE_LAT['hparam']
    hparam = hparam[np.argmax(reflat <= lmax):np.argmin(reflat >= lmin)]
    # NOTE: pretty sure anything over default valid_max is unphysical
    hparam[hparam > valid_max] = np.nan
    hparam = fill_invalid(hparam, method='nearest', inplace=True)
    return np.clip(_shrink(hparam, 'hparam', dpp), *model_bounds)


def load_hparam(
    dpp: float,
    model_bounds=(0.02, 0.13),
    valid_max=0.5,
    overwrite=False
) -> fits.HDUList:
    if INT_PATHS['hparam'].exists() and overwrite is False:
        return fits.open(INT_PATHS['hparam'])
    hparam = _manipulate_hparam(
        SRC_PATHS['hparam'], dpp, valid_max, model_bounds
    )
    hdul = _intermediate_hdul(hparam, 'hparam', dpp)
    hdul['HPARAM'].header['CLIPMIN'] = model_bounds[0]
    hdul['HPARAM'].header['CLIPMAX'] = model_bounds[1]
    hdul['HPARAM'].header['VALIDMAX'] = valid_max
    hdul.writeto(INT_PATHS['hparam'], overwrite=True)
    return hdul


def _hapke2albedo(albedo: np.ndarray) -> np.ndarray:
    """converts hapke parameter to albedo"""
    return -18 * albedo ** 2 + 3.39 * albedo + 0.00427


def _manipulate_albedo(sourcefile, dpp, valid_bounds, model_bounds):
    albedo = fits.open(sourcefile)[0].data.astype('f4')
    albedo[albedo <= valid_bounds[0]] = np.nan
    albedo[albedo >= valid_bounds[1]] = np.nan
    albedo = fill_invalid(albedo, inplace=True, method='nearest')
    albedo = np.clip(_shrink(albedo, 'albedo', dpp), *model_bounds)
    return _hapke2albedo(albedo)


def load_albedo(
    dpp: float,
    valid_bounds=(0, 0.189),
    model_bounds=(0.02, 0.16),
    overwrite: bool = False
) -> fits.HDUList:
    if INT_PATHS['albedo'].exists() and overwrite is False:
        return fits.open(INT_PATHS['albedo'])
    albedo = _manipulate_albedo(
        SRC_PATHS['albedo'], dpp, valid_bounds, model_bounds
    )
    hdul = _intermediate_hdul(albedo, 'albedo', dpp)
    for p, b in zip(('VALID', 'CLIP'), (valid_bounds, model_bounds)):
        hdul['ALBEDO'].header[f'{p}MIN'] = b[0]
        hdul['ALBEDO'].header[f'{p}MAX'] = b[1]
    hdul.writeto(INT_PATHS['albedo'], overwrite=True)
    return hdul


def _manipulate_dem(demtype, resolution):
    # These maps are centered at 180 longitude; everything
    # else is centered at 0; roll to match.
    dem = pdr.read(SRC_PATHS[demtype]).IMAGE.astype('f4')
    dem = np.roll(dem, dem.shape[1] // 2, axis=1)
    return _shrink(dem, demtype, resolution)


def load_dem(
    demtype: Literal["slope", "azimuth"], dpp: float, overwrite: bool = False
) -> fits.HDUList:
    if INT_PATHS[demtype].exists() and overwrite is False:
        return fits.open(INT_PATHS[demtype])
    hdul = _intermediate_hdul(_manipulate_dem(demtype, dpp), demtype, dpp)
    hdul.writeto(INT_PATHS[demtype], overwrite=True)
    return hdul


def make_filepaths(pipeline_data_path: Path, bundle_data_path: Path):
    srcmap_path = pipeline_data_path / "source_maps"
    intmap_path = pipeline_data_path / "intermediate_maps"
    tablepath = pipeline_data_path
    for k, v in SRCMAP_FNS.items():
        SRC_PATHS[k] = srcmap_path / v
        INT_PATHS[k] = intmap_path / f"{k}_processed.fits"
    SRC_PATHS["tb_table"] = tablepath / SRCTABLE_FN
    INT_PATHS["tb_modelmap_root"] = bundle_data_path
    notfound = [v.name for v in SRC_PATHS.values() if not v.exists()]
    if len(notfound) > 0:
        warnings.warn(
            f"Source map(s) not found: {', '.join(notfound)}. Map loading "
            f"will fail unless intermediate maps are present and "
            f"overwrite=False."
        )
    INT_PATHS["tb_modelmap_root"].mkdir(exist_ok=True)
    intmap_path.mkdir(exist_ok=True)


def load_source_maps(dpp: float, overwrite: bool) -> dict[str, fits.HDUList]:
    return {
        'hparam': load_hparam(dpp, overwrite=overwrite),
        'albedo': load_albedo(dpp, overwrite=overwrite),
        'slope': load_dem("slope", dpp, overwrite=overwrite),
        'azimuth': load_dem("azimuth", dpp, overwrite=overwrite)
    }


def prep_source_data(
    dpp: float, overwrite: bool = False,
) -> dict[str, np.ndarray]:
    source = load_source_maps(dpp, overwrite)
    out = {
        # note choice of hparam for lat/lon is somewhat arbitrary -- we just
        # need the smallest. it and albedo are the same; dem maps are bigger.
        'lon': source['hparam']['LONGITUDE'].data,
        'lat': source['hparam']['LATITUDE'].data,
        'hparam': source['hparam']['HPARAM'].data,
        'albedo': source['albedo']['ALBEDO'].data
    }
    # cut slope & azimuth to +/- 70 latitude to match hparam & albedo
    lat_ix = np.isin(source['slope']['LATITUDE'].data, out['lat'])
    for maptype in ('slope', 'azimuth'):
        out[maptype] = source[maptype][maptype.upper()].data[lat_ix]
    gmap(lambda h: h.close(), source.values())
    return out


def load_brightness_table(
    channel: Literal[1, 2, 3, 4], c_constant: float
) -> tuple[pa.Table, dict[str, float]]:
    tab = pq.read_table(
        SRC_PATHS['tb_table'],
        filters=[
            ('freq', '=', np.float32(CHANNEL_FREQ[channel])),
            ('c', '=', c_constant)
        ],
        columns=['lat', 'albedo', 'h', 'ltst', 'tb']
    )
    return tab, pq_metadict(SRC_PATHS['tb_table']) | {'c': c_constant}


def dem_ltst_shift(
    slope: np.ndarray, azimuth: np.ndarray, max_shift: float = 1.5
) -> np.ndarray:
    return np.clip(
        slope * np.sin(np.radians(azimuth)) * 24 / 360, -max_shift, max_shift
    )


def _find_timebounds(high, low, ltst_precision, offset):
    if low < 0 and high > offset:
        raise ValueError("Something has gone wrong with time correction.")
    elif low < 0:
        bounds = (round(offset + low, ltst_precision), high)
        wrapped = True
    elif high > offset:
        bounds = (low, round(high - offset, ltst_precision))
        wrapped = True
    else:
        bounds, wrapped = (low, high), False
    op = pc.or_ if wrapped is True else pc.and_
    return bounds, op, wrapped


def prep_timebin_table(
    tab: pa.Table,
    source: dict[str, np.ndarray],
    timebin: tuple[float, float],
    max_shift: float,
    ltst_precision: int
) -> tuple[pa.Table, np.ndarray]:
    ltst = np.mean(timebin) + dem_ltst_shift(
        source['slope'], source['azimuth'], max_shift
    )
    low, high = ltst.min(), ltst.max()
    offset = round(24 - 1 / (10 * ltst_precision), ltst_precision)
    bounds, op, wrapped = _find_timebounds(high, low, ltst_precision, offset)
    matches = op(
        pc.greater_equal(tab['ltst'], bounds[0]),
        pc.less_equal(tab['ltst'], bounds[1])
    )
    tab = tab.filter(matches)
    if low < 0:
        ltstcol = pc.choose(
            pc.greater_equal(tab['ltst'], low + 24),
            tab['ltst'],
            pc.subtract(tab['ltst'], pa.scalar(24, pa.float32())),
        )
    elif high > offset:
        ltstcol = pc.choose(
            pc.less_equal(tab['ltst'], high - offset),
            tab['ltst'],
            pc.add(tab['ltst'], pa.scalar(24, pa.float32())),
        )
    if wrapped is True:
        # noinspection PyUnboundLocalVariable
        tab = tab.set_column(tab.column_names.index('ltst'), 'ltst', ltstcol)
    return tab, ltst


def run_timebin(
    tb_tab: pa.Table,
    ltst: np.ndarray,
    lat: np.ndarray,
    source: dict[str, np.ndarray],
):
    variables = {
        'albedo': source['albedo'],
        'h': source['hparam'],
        # note that model cares only about absolute latitude
        'lat': abs(lat),
        'ltst': ltst
    }
    dependent_variable = 'tb'
    tb_tab = tb_tab.sort_by([(v, "ascending") for v in variables.keys()])
    # note that this only works because we sorted it in the correct order
    # first! That sort basically turned it into a last-index-fastest raster.
    varpts = [
        pc.unique(tb_tab[c]).to_numpy()
        for c in variables.keys()
    ]
    depgrid = tb_tab[dependent_variable].to_numpy().reshape(
        [varpts[i].size for i in range(len(varpts))]
    )
    samplegrid = np.dstack(
        [
            # note that we nearest-neighbor interpolated over bad values in the
            # source products while constructing the intermediate products.
            # We also clipped them, and shouldn't actually need to clip again
            # here; this is an abundance of caution.
            np.clip(v, varpts[i].min(), varpts[i].max())
            for i, v in enumerate(variables.values())
        ]
    )
    # NOTE: is there a way to prevent this from internally casting to float64?
    return interpn(varpts, depgrid, samplegrid).astype('f4')


def run_channel(
    channel: Literal[1, 2, 3, 4],
    timebins: Sequence[tuple[float, float]],
    resolution: float = 1 / 32,
    overwrite_intermediate_maps: bool = False,
    max_time_correction: float = 1.5,
    ltst_precision: int = 1,
    c_constant: float = 2.5
):
    print(f"making interpolated TB model map for channel {channel}...")
    # create or load all source-data maps as required. these serve to provide
    # georegistered interpolation points within the parameter space table
    # loaded by load_brightness_table().
    source = prep_source_data(resolution, overwrite_intermediate_maps)
    # construct lat array. note that the model only cares about absolute
    # latitude, and not at all about longitude. choice of the 'hparam' map is
    # arbitrary -- they should all be the same shape.
    # TODO: can cut the 'f4' later
    lat = np.tile(source['lat'], (source['hparam'].shape[1], 1)).T.astype('f4')
    # load modeled brightness temperature table
    # TODO, maybe: c_value will be variable if we can do loss tangent stuff.
    # TODO: we may also use the slope/az lat correction at some point.
    print("loading table...")
    tab, meta = load_brightness_table(channel, c_constant)
    meta |= {
        'PPD': int(1 / resolution),
        'CHANNEL': channel,
        'MAXTCORR': max_time_correction
    }
    hdus = [fits.PrimaryHDU()]
    for timebin in timebins:
        print(f"mapping channel {channel} Tb model in timebin {timebin}...")
        tb_tab, timebin_ltst = prep_timebin_table(
            tab, source, timebin, max_time_correction, ltst_precision
        )
        hdu = to_scaled_hdu(
            run_timebin(tb_tab, timebin_ltst, lat, source),
            # NOTE: as elsewhere, if we do this with noninteger timebins, this
            #  will need to change
            name=f"TBMOD_{'_'.join(map(str, timebin))}"
        )
        hdu.header = fits.Header(dict(hdu.header) | meta)
        hdu.header['PTYPE'] = 'TBMOD'
        hdus.append(hdu)
        del tb_tab
    del tab
    hdus.append(fits.ImageHDU(lat[:, 0], name="LATITUDE"))
    hdus.append(fits.ImageHDU(source['lon'].astype('f4'), name="LONGITUDE"))
    print(f"writing FITS file for channel {channel}...")
    ppdstr = f"{int(1 / resolution)}ppd"
    fits.HDUList(hdus).writeto(
        INT_PATHS["tb_modelmap_root"] / f"t{channel}_tbmod_{ppdstr}.fits",
        overwrite=True
    )
    del hdus
    print(f"done with channel {channel}.\n")
