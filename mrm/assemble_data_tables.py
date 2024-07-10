"""
This script further processes intermediate concatenated versions of the source
ASCII tables. It produces an intermediate parquet file used by run_deconv.py
as well as FITS tables intended for inclusion in the bundle's data collection.

Various aspects of its behavior can be controlled by modifying constants
described in the body of the script.

Note that this script does not generate PDS4 labels for the FITS table files,
which must be modified manually if you wish to create valid products of
different lengths etc.
"""
from pathlib import Path

from astropy.io import fits
from astropy.table import Table
import numpy as np
import pandas as pd
import pyarrow as pa
from mrm.processing.tables import (
    channel_offset,
    ESETTINGS,
    et2lst_map,
    jd_tdb_map,
    latmodel_outliers,
    latcut_outliers,
    latsrf_map,
    low_temp,
    bad_orbits,
    exclude_orbits,
    PROFILER,
    spkcpo_map,
    small_orbits,
    timedupes,
    utc2et_map,
)
from lhorizon import LHorizon
from lhorizon.constants import LUNAR_RADIUS
from lhorizon.lhorizon_utils import hats
from pyarrow import parquet as pq, compute as pac


# *****CONFIGURATION*****
# which orbiters should we generate tables for? if you haven't already made
# tables for each orbiter, this should always be ('ce1', 'ce2')
ORBITER_NAMES = ('ce1', 'ce2')
# which exclusions (defined in ephemeris_assembly.py) shall we perform?
EXCLUSIONS = (
    low_temp,
    channel_offset,
    small_orbits,
    timedupes,
    bad_orbits,
    exclude_orbits,
    latmodel_outliers,
    latcut_outliers
)
# the exclusion functions defined in ephemeris_assembly.py
# reference values inserted here into ESETTINGS.
ESETTINGS['latmodel_cutoff'] = 8
ESETTINGS['min_orbit_size'] = 100
ESETTINGS['min_temp'] = 34
ESETTINGS['max_channel_offset'] = 75
ESETTINGS['tile_vote_centile'] = 98
# number of threads for CPU-bound operations (None disables multithreading)
THREADS = 8
# kwargs used for spot checks against Horizons
HKWARGS = {
    'query_type': 'VECTORS',
    'query_options': {'refsystem': 'ICRF', 'vec_corr': 'LT+S'}
}
# where is the metakernel file in which we have defined the kernels
# SPICE should use for calculations?
ESETTINGS["metakernel"] = "lhorizon_metakernel.tm"
# where are the parquet files generated from the L2C tables?
INPATH = Path(__file__).absolute().parent / "pipeline_data/mrm_l2c_tables"
# where shall we write intermediate split tables?
SPLITPATH = (
    Path(__file__).absolute().parent.parent / "pipeline_data/split_tables"
)
# where shall we write the intermediate concatenated table?
CATPATH = Path(__file__).absolute().parent / "pipeline_data/concatenated_table"
# where shall we write the FITS tables for archival?
PDSPATH = Path(__file__).absolute().parent / "bundle/data/"
# which orbits / orbit ranges do we want to call 'bad' or exclude for
# calibration reasons?
ORBIT_EXCLUSIONS = {
    'ce1': {
        'bad_orbits': np.array([243, 4856]),
        'exclude_orbits': np.array([], dtype=int)
    },
    'ce2': {
        'bad_orbits': np.array(
            [102,103,112,113,114,169,170,171,172,173,174,275,276,277,278,279,
             280,281,282,283,284,285,286,287,288,289,290,291,780,781,782,783,
             792,793,794,856,857,858,859,1388,1389,1390,1502,1839,2314,2315,
             2316,2317,2318,2319,2320,2321,2642,2643,2644,2645,2646,2647,2648,
             2649,2650,2651,2652,2653,2654]),
        'exclude_orbits': np.arange(1470, 1640)
    }
}
"""Fields of intermediate parquet table."""
PA_FIELDS = (
    ('orbit', pa.uint16()),
    ('utc', pa.timestamp('ms')),
    ('jd', pa.float64()),
    ('tdb', pa.string()),
    ('et', pa.float64()),
    ('ltst', pa.float32()),
    *((f't{i}', pa.float32()) for i in range(1, 5)),
    ('theta', pa.float32()),
    ('phi', pa.float32()),
    ('lat', pa.float32()),
    ('lon', pa.float32()),
    ('d', pa.float32()),
    *((f'b{ax}', pa.float64()) for ax in ('x', 'y', 'z')),
    *((f'c{ax}', pa.float64()) for ax in ('x', 'y', 'z')),
    *((f'm{ax}', pa.float64()) for ax in ('x', 'y', 'z')),
    *((f'n{ax}', pa.float64()) for ax in ('x', 'y', 'z')),
    *((f'o{ax}', pa.float32()) for ax in ('x', 'y')),
    ('flag', pa.uint16())
)
PA_SCHEMA = pa.schema(PA_FIELDS)


def apply_mission_flags(orbiter, orbname):
    mission_ok = "0" if orbname == 'ce2' else "0X000000"
    mission_flagged = pac.invert(
        pac.equal(orbiter['flag'], mission_ok)
    ).to_numpy()
    flagsum = mission_flagged.sum()
    print(
        f"{flagsum} samples flagged by mission "
        f"({round(flagsum / len(orbiter) * 100, 2)}%) "
        f"(flag value 1)"
    )
    return mission_flagged


def process_orbiter(orbname):
    print(f"*****processing {orbname}*****")
    PROFILER.reset()
    ESETTINGS["bad_orbits"] = ORBIT_EXCLUSIONS[orbname]["bad_orbits"]
    ESETTINGS["exclude_orbits"] = ORBIT_EXCLUSIONS[orbname]["exclude_orbits"]
    # load table created from L2C data.
    # note that PROFILER is simply a performance profiler.
    print("loading table...")
    with PROFILER.context("pa_io"):
        input_table_path = INPATH / f'{orbname}_l2c_complete.parquet'
        # always drop spatially invalid data not captured by flagging
        orbtab = pq.read_table(
            input_table_path,
            filters=[("lon", "<=", 360), ("lat", "<=", 90), ("d", "<=", 999.9)]
        )
    # flag 'bad' data
    print("applying quality flags...")

    with PROFILER.context("exclusions"):
        ESETTINGS['flagarr'] = np.zeros(len(orbtab), 'u2')
        # always flag rows defined as bad by the original data providers
        mission_flagged = apply_mission_flags(orbtab, orbname)
        ESETTINGS['flagarr'][mission_flagged] = 1
        for e in EXCLUSIONS:
            flagged, flagval = e(orbtab)
            ESETTINGS['flagarr'][flagged] += flagval
            flagsum = len(pac.indices_nonzero(flagged))
            print(
                f"{flagsum} samples flagged as {e.__name__} "
                f"({round(flagsum / len(orbtab) * 100, 2)}%) "
                f"(flag value {flagval})"
            )
        flag_total = len(np.nonzero(ESETTINGS['flagarr'] != 0)[0])
        print(
            f"{flag_total} total samples flagged "
            f"({round(flag_total / len(orbtab) * 100, 2)}%)"
        )

        print("rebuilding table...")
        with PROFILER.context("array_shaping"):
            orbtab = orbtab.to_pandas()
            orbtab['flag'] = ESETTINGS['flagarr']
            # convert UTC times to strings for use by astropy and SPICE;
            # also pre-chunk them to map across processes in later steps
            nchunks = 1 if THREADS is None else THREADS
            utc_chunks = np.array_split(
                orbtab['time'].astype(str).values, nchunks
            )
        # compute ET (time scale used by SPICE) for each sample.
        # as above, although offsets from vectorized methods
        # are probably below our resolution, let's introduce as little slop
        # into our calculations whatsoever by using SPICE's utc2et.
        print("computing Ephemeris Time...")
        with PROFILER.context("SPICE time conversion"):
            et_time = utc2et_map(utc_chunks, THREADS)
        # recompute orbiter TDB and JD values with astropy
        # to ensure millisecond-level precision
        print("computing JD and TDB...")
        with PROFILER.context('astropy time conversion'):
            jd, tdb = jd_tdb_map(utc_chunks, THREADS)

        # compute local solar time for each sample using SPICE
        print("computing local time...")
        with PROFILER.context('ltst computation'):
            ltst = et2lst_map(et_time, orbtab['lon'], THREADS)
        # find vectors for orbiter boresight & orbiter
        # in MOON_ME (body-fixed) frame
        # at least using the ELLIPSOID method, this is essentially
        # no different from a naive spherical-to-cartesian method,
        # but dotting the is and crossing the ts, etc.
        print("finding selenocenter-boresight vectors...")
        with PROFILER.context('latsrf'):
            selenocenter_boresight_vectors = latsrf_map(
                et_time, orbtab['lon'], orbtab['lat'], THREADS
            )
        # simple geometry to get orbiter position in body-fixed frame.
        # assuming the Moon is a smooth sphere and orbiter is always looking
        # perfectly nadir, the position vector from selenocenter to orbiter
        # is always in the direction of the surface normal vector.
        print("projecting body-fixed vectors...")
        with PROFILER.context("projection"):
            selenocenter_orbiter_vectors = np.einsum(
                'ij,i->ij',
                hats(selenocenter_boresight_vectors),
                orbtab['d'] + LUNAR_RADIUS / 1000
            )
            surface_orbiter_vectors = (
                selenocenter_boresight_vectors - selenocenter_orbiter_vectors
            )
            # quick absurdity check
            assert np.allclose(
                np.linalg.norm(surface_orbiter_vectors, axis=1),
                orbtab['d']
            )
        # transform coordinates: MOON_ME (body-fixed) -> inertial (ECLIPJ2000)
        print("computing inertial coordinates...")
        with PROFILER.context('SPICE position computation'):
            bore_sun_vectors = spkcpo_map(
                et_time, selenocenter_boresight_vectors, THREADS
            )
            print('boresight vectors computed')
            orbiter_sun_vectors = spkcpo_map(
                et_time, selenocenter_orbiter_vectors, THREADS
            )
            print('orbiter vectors computed')
            moon_sun_vectors = spkcpo_map(et_time, [0, 0, 0], THREADS)
            print('moon vectors computed')
        # spot-check against horizons
        print("spot-checking against JPL Horizons...")
        indices = np.random.randint(0, len(et_time), 5)
        for ix in indices:
            horizons_location = {
                'lat': orbtab['lat'].iloc[ix],
                'lon': orbtab['lon'].iloc[ix],
                'elevation': 0,
                'body': 301
            }
            check = LHorizon(
                target='10', origin=horizons_location, epochs=tdb[ix],
                **HKWARGS
            )
            horizons_position = check.dataframe()[['X', 'Y', 'Z']].values
            offset = horizons_position - bore_sun_vectors[ix]
            print(np.linalg.norm(offset))
            # consider less than 50m offset fine -- average is more like 20m,
            # and probably mostly due to floating-point error
            assert np.linalg.norm(offset) < 0.05, "offset too large"
        print("concatenating table segments...")
        segments = [
            orbtab,
            pd.DataFrame(bore_sun_vectors, columns=['bx', 'by', 'bz']),
            pd.DataFrame(orbiter_sun_vectors, columns=['cx', 'cy', 'cz']),
            pd.DataFrame(moon_sun_vectors, columns=['mx', 'my', 'mz']),
            pd.DataFrame(selenocenter_boresight_vectors,
                         columns=['nx', 'ny', 'nz']),
        ]
        # concatenate segments, then add orthographic projection, ET, TDB,
        # and LTST to ephemeris
        full_eph = pd.concat(segments, axis=1).rename(columns={'time': 'utc'})
        full_eph['et'] = et_time
        full_eph['ltst'] = ltst
        full_eph['tdb'] = tdb
        print("building pyarrow table...")
        with PROFILER.context("array_shaping"):
            full_eph = full_eph.reindex(columns=PA_SCHEMA.names)
            table = pa.Table.from_pandas(full_eph, preserve_index=False)
            del full_eph
            table = table.cast(pa.schema(PA_SCHEMA))
        print("writing parquet table...")
        SPLITPATH.mkdir(exist_ok=True, parents=True)
        with PROFILER.context('pa_io'):
            pq.write_table(
                table, SPLITPATH / f"{orbname}_ephemeris.parquet"
            )
        print(PROFILER)


def concatenate_split_tables():
    CATPATH.mkdir(exist_ok=True, parents=True)
    tabs = []
    for o in ORBITER_NAMES:
        tab = pq.read_table(SPLITPATH / f"{o}_ephemeris.parquet")
        tab = tab.append_column('orbiter', pa.array(np.full(len(tab), o)))
        tabs.append(tab)
    pq.write_table(
        pa.concat_tables(tabs), CATPATH / 'ce_mrm_concatenated.parquet'
    )


def colix(tab, col):
    return tab.column_names.index(col)


def setcol(tab, col, arr):
    return tab.set_column(colix(tab, col), col, arr)


def write_pds_fits_table(orbiter, input_table):
    tab = pq.read_table(input_table, filters=[('orbiter', '=', orbiter)])
    tab = setcol(tab, "utc", tab['utc'].cast(pa.string()))
    for ax in ("x", "y", "z"):
        tab = setcol(tab, f"n{ax}", tab[f"n{ax}"].cast(pa.float32()))
    tab = tab.drop(["jd", "tdb", "ox", "oy", "theta", "phi", "orbiter"])
    t1 = Table()
    for col in tab.column_names:
        print(f"{col}...", end="")
        colarr = tab[col].to_numpy()
        if colarr.dtype.char == 'O':
            colarr = colarr.astype(str)
        t1[col] = colarr
    hdu = fits.BinTableHDU(t1, name="table")
    hdul = fits.HDUList([fits.PrimaryHDU(), hdu])
    Path("products/pds/temp").mkdir(exist_ok=True, parents=True)
    hdul.writeto(PDSPATH / f'{orbiter}_mrm.fits', overwrite=True)


def write_fits_tables():
    input_table = CATPATH / 'ce_mrm_concatenated.parquet'
    PDSPATH.mkdir(exist_ok=True, parents=True)
    for o in ORBITER_NAMES:
        write_pds_fits_table(o, input_table)


def assemble_tables():
    for o in ORBITER_NAMES:
        process_orbiter(o)
    print("concatenating split tables...")
    concatenate_split_tables()
    print("writing FITS archival tables...")
    write_fits_tables()


if __name__ == "__main__":
    assemble_tables()
