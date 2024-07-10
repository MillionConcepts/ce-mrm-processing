"""
This handler script produces data (.fits) and browse (.png) files for the
map products in the data collection and the browse products in the browse
collection. It requires the "base" deconvolved maps produced by run_deconv.py
along with the model outputs produced by running execute_heat1d.py,
assemble_heat1d_outputs.py, assemble_brightness_table.py, and
map_brightness.py in order. Various behaviors can be modified by editing
constants in the body of the script or passing parameters to the command-line
execution.
"""

import re
from numbers import Number
from pathlib import Path
from typing import Literal, Sequence

from astropy.io import fits
from cytoolz import dissoc, get
from dustgoggles.dynamic import exc_report
from dustgoggles.structures import listify, MaybePool
import fire
from hostess.monitors import Stopwatch
from marslab.imgops.imgutils import centile_clip, std_clip
from marslab.imgops.pltutils import (
    attach_axis,
    dpi_from_image,
    set_colorbar_font,
)
import matplotlib as mpl
import matplotlib.font_manager as mplf
import matplotlib.pyplot as plt
from more_itertools import all_equal
import numpy as np
import pandas as pd
from mrm.processing.projection import coef_det
from quickseries.simplefit import fit
from rich import print as rprint

from mrm.processing.hdus import to_scaled_hdu, decimal_to_hourstring
from mrm.processing.tables import lat_tb_model
from mrm.shared.style import NoShow, register_cmaps
from mrm.shared.constants import CHANNEL_FREQ
from mrm.shared.style import EASYDARK, FONT_PATH

DIVERGING_CMAPS = ("orange_teal", "neon")

register_cmaps()


def load_and_cut(inpath, orb, channel, max_abslat, variable_meta, ppd):
    print("loading...")
    hdulists = []
    for f in sorted(filter(lambda p: p.suffix == ".fits", inpath.iterdir())):
        if (hdul := fits.open(f))["PRIMARY"].header["CHANNEL"] != channel:
            continue
        if hdul["PRIMARY"].header["ORBITER"].strip().lower() != orb.lower():
            continue
        if hdul["PRIMARY"].header["PPD"] != ppd:
            continue
        hdulists.append(hdul)
    hdulists = sorted(
        hdulists, key=lambda h: h[0].header["OLTSTMIN"], reverse=True
    )
    out = {"headers": [h[0].header for h in hdulists]}
    # ensure input files all came from pipeline runs with compatible settings,
    # for configurable values of 'compatible'
    pin_equal = tuple(
        map(lambda h: dissoc(dict(h), *variable_meta), out["headers"])
    )
    if not all_equal(pin_equal):
        raise ValueError("Files appear to be from mismatched runs.")
    # longitude bins should all be equal -- if not, it means that there was
    # much worse (or differently localized) coverage than is expected in any
    # of these maps. This might be a bad assumption if we start making quite
    # small local maps using the same method, but we probably shouldn't anyway.
    lonbins = np.vstack([h["LONBINS"].data for h in hdulists])
    if not np.unique(np.diff(lonbins, axis=0))[0] == 0:
        raise ValueError("Misaligned lonbins")
    # latitude bins will not be stackable initially (because they don't wrap,
    # so small differences in availability can matter)
    valid_latbins = [
        np.nonzero(abs(h["LATBINS"].data) <= max_abslat)[0] for h in hdulists
    ]
    try:
        valid_latbins = np.vstack(valid_latbins)
    except ValueError:
        raise ValueError("Missing latbins")
    latbins = np.vstack(
        [h["LATBINS"].data[v] for h, v in zip(hdulists, valid_latbins)]
    )
    if not np.unique(np.diff(latbins, axis=0))[0] == 0:
        raise ValueError("Misaligned latbins")
    out["val_arrays"] = [
        h[1].data[vl, :] for h, vl in zip(hdulists, valid_latbins)
    ]
    out["std_arrays"] = [
        h[2].data[vl, :] for h, vl in zip(hdulists, valid_latbins)
    ]
    return out | {"latbins": latbins[0], "lonbins": lonbins[0]}


def assemble_base_hdul(headers, lonbins, latbins, val_arrays, std_arrays):
    """Build single-channel stacked-time TEMP/STD HDUList from loaded data"""
    val_hdus, std_hdus = [], []
    ch = headers[0]['CHANNEL']
    for i in range(len(headers)):
        h, v, s = map(list.pop, (headers, val_arrays, std_arrays))
        tstr = decimal_to_hourstring(h["OLTSTMIN"], h["OLTSTMAX"])
        print(f"processing T{ch} TEMP/STDEV {tstr}...", end="", flush=True)
        val_hd = fits.Header(dict(h) | {"PTYPE": "TEMP", "BUNIT": "K"})
        val_hdus.append(to_scaled_hdu(v, name=f"TEMP_{tstr}", header=val_hd))
        std_hd = fits.Header(dict(h) | {"PTYPE": "STDEV", "BUNIT": "K"})
        std_hdus.append(to_scaled_hdu(s, name=f"STDEV_{tstr}", header=std_hd))
    # NOTE: for an equirectangular projection, creating full backplanes from
    #  lat/lon is silly, so we're leaving them as 1D for space.
    ax_hdus = [
        fits.ImageHDU(latbins, name="LATITUDE"),
        fits.ImageHDU(lonbins, name="LONGITUDE"),
    ]
    for hdu in ax_hdus:
        hdu.header["BUNIT"] = "deg"
    return fits.HDUList([fits.PrimaryHDU(), *val_hdus, *std_hdus, *ax_hdus])


def write_base_hdul(
    slice_path, map_path, orbiter, channel, max_abslat, variable_meta, ppd
):
    """
    Create and write single-channel stacked-time TEMP/STD FITS file from
    single channel/timebin FITS files produced by deconvolution pipeline
    """
    data = load_and_cut(
        slice_path, orbiter, channel, max_abslat, variable_meta, ppd
    )
    hdul = assemble_base_hdul(**data)
    orb, ppd = (hdul[1].header[k] for k in ("ORBITER", "PPD"))
    path = map_path / f"{orb.lower()}_t{channel}_temp_{ppd}ppd.fits"
    print(f"writing to {path}...", end="", flush=True)
    hdul.writeto(path, overwrite=True)


def best_tick_positions(
    arr: np.ndarray,
    targets: Sequence[Number],
    tolerance: float = 0,
) -> list[np.ndarray]:
    positions = []
    for t in targets:
        valid = np.nonzero(abs(arr - t) <= tolerance)[0]
        if valid.size == 0:
            raise ValueError(f"No match within {tolerance} for {t}.")
        positions.append(valid[valid.size // 2])
    return positions


def _load_raveled(fpath, ix) -> tuple[fits.Header, np.ndarray]:
    hdu: fits.ImageHDU = fits.open(fpath)[ix]
    return hdu.header, hdu.data.ravel()


def _simple_lat_fit(fitlat, n_fit_points, temp):
    fit_df = (
        pd.DataFrame({"lat": fitlat, "temp": temp})
        .dropna(axis=0)
        .sample(n_fit_points)
        .reset_index(drop=True)
    )
    fit_params, _ = fit(
        lat_tb_model,
        [fit_df["lat"]],
        fit_df["temp"],
        guess=[231, 0.35]
    )
    det = coef_det(lat_tb_model(fit_df['lat'], *fit_params), fit_df['temp'])
    return fit_params, det


def assemble_latshift_hdul(inpath, orb, channel, ppd, n_fit_points=int(2e6)):
    """
    Load data from per-channel stacked-time TEMP/STDEV FITS files and assemble
    it into a LATSHIFT HDUList.

    n_fit_points is size of random per-timebin subsample for latitudinal
    fitting. The fit can otherwise take 20+ seconds per slice at 32ppd,
    which is annoying, and a random subsample this size should be fine.
    """
    fpath = inpath / f"{orb.lower()}_t{channel}_temp_{ppd}ppd.fits"
    temp_hdu_indices = [
        r[0] for r in fits.open(fpath).info([]) if r[1].startswith("TEMP")
    ]
    shape = get(["NAXIS2", "NAXIS1"], fits.open(fpath)[1].header)
    # latitude bins
    baselat = np.radians(fits.open(fpath)["LATITUDE"].data).astype("f4")
    # latitude bins indexed to match raveled temps in fit process
    fitlat = baselat[np.indices(shape)[0].ravel()]
    lathdus = []
    for ix in temp_hdu_indices:
        header, temp = _load_raveled(fpath, ix)
        extname = header["EXTNAME"].replace("TEMP", "LATSHIFT")
        header["PTYPE"] = "LATSHIFT"
        print(f"processing T{channel} {extname}...", end="", flush=True)
        fit_params, det = _simple_lat_fit(fitlat, n_fit_points, temp)
        header["MODELDET"] = np.float32(det)
        header["MPARAM_A"], header["MPARAM_B"] = fit_params
        difference = (temp - lat_tb_model(fitlat, *fit_params)).reshape(shape)
        del temp
        lathdus.append(to_scaled_hdu(difference, name=extname, header=header))
    return fits.HDUList(
        [
            fits.PrimaryHDU(),
            *lathdus,
            fits.open(fpath)["LATITUDE"],
            fits.open(fpath)["LONGITUDE"],
        ]
    )


def write_latshift_hdul(map_path, orbiter, channel, ppd):
    """Build and write single-channel stacked-time LATSHIFT FITS file"""
    hdul = assemble_latshift_hdul(map_path, orbiter, channel, ppd)
    path = map_path / f"{orbiter.lower()}_t{channel}_latshift_{ppd}ppd.fits"
    print(f"writing to {path}...", end="", flush=True)
    hdul.writeto(path, overwrite=True)


def assemble_datminus_hdul(map_path, orb, channel, ppd):
    """
    Load data from per-channel stacked-time TEMP/STDEV and TBMOD FITS files,
    difference them, and assemble them into a DATMINUS HDUL.
    """
    obs = map_path / f"{orb.lower()}_t{channel}_temp_{ppd}ppd.fits"
    mod = map_path / f"t{channel}_tbmod_{ppd}ppd.fits"
    assert obs.exists(), f"missing TEMP for t{channel} / {ppd}ppd / ce{orb}"
    assert mod.exists(), f"missing TBMOD for t{channel} / {ppd}ppd"
    obslon, modlon = (fits.open(f)['LONGITUDE'].data for f in (obs, mod))
    assert np.all(obslon == modlon), "mismatched longitude bins."
    # Find latitude bin indices for temp arrays that match model lat.
    # NOTE: this assumes that we rendered temp data to at least model latitude
    #  range. In current default settings this is true: +/- 75 lat vs. +/- 70
    #  lat, global longitude. If this is not true at some point in the future,
    #  this will need to change.
    lat_ix = np.isin(*[fits.open(f)['LATITUDE'].data for f in (obs, mod)])
    temp_hdu_indices = [
        r[0] for r in fits.open(obs).info([]) if r[1].startswith("TEMP")
    ]
    mod_hdu_indices = [
        r[0] for r in fits.open(mod).info([]) if r[1].startswith("TBMOD")
    ]
    dathdus = []
    for tempix, modix in zip(temp_hdu_indices, mod_hdu_indices):
        temphdu = fits.open(obs)[tempix]
        header, temp = temphdu.header, temphdu.data[lat_ix]
        del temphdu
        extname = header['EXTNAME'].replace('TEMP', 'DATMINUS')
        print(f"processing T{channel} {extname}...", end="", flush=True)
        header['PTYPE'] = 'DATMINUS'
        dathdu = to_scaled_hdu(
            temp - fits.open(mod)[modix].data, header=header, name=extname
        )
        dathdus.append(dathdu)
    return fits.HDUList(
        [
            fits.PrimaryHDU(),
            *dathdus,
            fits.open(mod)["LATITUDE"],
            fits.open(mod)["LONGITUDE"],
        ]
    )


def write_datminus_hdul(map_path, orbiter, channel, ppd):
    """Build and write single-channel stacked-time DATMINUS FITS file"""
    hdul = assemble_datminus_hdul(map_path, orbiter, channel, ppd)
    path = map_path / f"{orbiter.lower()}_t{channel}_datminus_{ppd}ppd.fits"
    print(f"writing to {path}...", end="", flush=True)
    hdul.writeto(path, overwrite=True)


def map_to_browse(
    arr: np.ndarray,
    latbins: np.ndarray,
    lonbins: np.ndarray,
    outpath: Path,
    latgrid: np.ndarray = np.arange(-60, 90, 30),
    longrid: np.ndarray = np.arange(-150, 180, 50),
    dpi_downscale_factor: int = 3,
    fontproperties=mplf.FontProperties(),
    title=None,
    clip=None,
    cmap=plt.get_cmap("Greys_r"),
    normtype="standard",
    **imshow_kwargs,
):
    with NoShow(), plt.style.context(EASYDARK):
        fig, ax = plt.subplots()
        clip = {} if clip is None else clip
        if clip.get("type") == "centile":
            arr = centile_clip(arr, clip["value"])
        elif clip.get("type") == "sigma":
            arr = std_clip(arr, clip["value"])
        elif clip.get("type") == "absolute":
            arr = np.clip(arr, *clip["value"])
        elif clip.get("type") is not None:
            raise ValueError(f"Unknown clip type {clip['type']}")
        if normtype == "centered" and clip.get("type") == "absolute":
            norm = mpl.colors.TwoSlopeNorm(0, *clip["value"])
        elif (
            normtype == "centered"
            and cmap.name.replace("_r", "") in DIVERGING_CMAPS
        ):
            arr = arr - np.nanmean(arr)
            norm = mpl.colors.TwoSlopeNorm(0, np.nanmin(arr), np.nanmax(arr))
        elif clip.get("type") == "absolute":
            norm = mpl.colors.Normalize(
                vmin=clip['value'][0], vmax=clip['value'][1]
            )
        else:
            norm = mpl.colors.Normalize()
        mappable = ax.imshow(arr, norm=norm, cmap=cmap, **imshow_kwargs)
        ax.set_xticks(
            best_tick_positions(lonbins, longrid),
            longrid,
            fontproperties=fontproperties,
        )
        ax.set_yticks(
            best_tick_positions(latbins, latgrid),
            latgrid,
            fontproperties=fontproperties,
        )
        ax.grid(color=(1, 1, 1), lw=1, alpha=0.35)
        colorbar = plt.colorbar(mappable, attach_axis(ax, size="3%"))
        set_colorbar_font(colorbar, fontproperties)
        if title is not None:
            title_fp = fontproperties.copy()
            title_fp.set_size(fontproperties.get_size() * 1.5)
            ax.set_title(title, fontproperties=title_fp)
        fig.tight_layout(pad=0)
        # return fig
        fig.savefig(
            outpath,
            dpi=dpi_from_image(fig) // dpi_downscale_factor,
            pad_inches=0.05,
            bbox_inches="tight",
        )
        plt.close(fig)


def make_browse_products(
    *,
    channel,
    map_path,
    browse_path,
    ptypes,
    cmaps,
    norms,
    clips,
    fontproperties,
    dpi_downscale,
    orbiter
):
    for path in filter(
        lambda p: re.search(rf"t{channel}.*\.fits", p.name),
        map_path.iterdir(),
    ):
        if path.name.startswith("ce") and not path.name.startswith(orbiter):
            continue
        latbins = fits.open(path)["LATITUDE"].data
        lonbins = fits.open(path)["LONGITUDE"].data
        for rec in filter(
            lambda r: re.search("|".join(ptypes), r[1]),
            fits.open(path).info([]),
        ):
            print(f"making browse image(s) for {path.name} {rec[1]}")
            header = fits.open(path)[rec[0]].header
            orb, channel, extname = (
                header.get(k) for k in ("ORBITER", "CHANNEL", "EXTNAME")
            )
            if (ptype := header.get("PTYPE")) == "TEMP":
                title = f"{orb} Tb (K) "
            elif ptype == "STDEV":
                title = (
                    f"Biased estimate of standard deviation of {orb} Tb (K) "
                )
            elif ptype == "LATSHIFT":
                title = f"{orb} Tb - 1st-order latitudinal trend (K) "
            elif ptype == "DATMINUS":
                title = f"{orb} Observed Tb - modeled Tb (K) "
            elif ptype == "TBMOD":
                title = f"modeled Tb (K) "
            else:
                raise NotImplementedError(f"unknown ptype {ptype}")
            title += f", MRM channel {channel} ({CHANNEL_FREQ[channel]} GHz) "
            numbers = re.search(rf"_(\d+)_(\d+)", extname).groups()
            title += (
                f", {numbers[0].zfill(2)}:00 to "
                f"{numbers[1].zfill(2)}:00 LTST. "
            )
            if (clip := clips.get(ptype, {})).get("type") == "centile":
                title += (
                    f"Clipped at centiles "
                    f"{'/'.join(map(str, clip['value']))}. "
                )
            elif clip.get("type") == "sigma":
                title += f"Clipped at +/- {clip['value']} sigma."
            elif clip.get("type") == "absolute":
                title += f"Clipped at {'/'.join(map(str, clip['value']))}."
            orbstr = "" if orb is None else f"{orb.upper()}_"
            for cn in (cmap_names := listify(cmaps[ptype])):
                if len(cmap_names) == 1:
                    fn = f"{orbstr}t{channel}_{extname}.png"
                else:
                    fn = (
                        f"{orbstr}t{channel}_{extname}_"
                        f"{cn.replace('_r', '')}.png"
                    )
                if (
                    cn.replace("_r", "") in DIVERGING_CMAPS
                    and clip.get("type") != "absolute"
                ):
                    suffix = " Zeroed at mean."
                else:
                    suffix = ""
                cmap = plt.colormaps.get_cmap(cn)
                cmap.set_bad((0, 0, 0), alpha=0)
                folder = browse_path / ptype.lower()
                if orb is not None:
                    folder = folder / orb
                folder = folder / f"t{channel}"
                if len(cmap_names) > 1:
                    folder = folder / cn.replace("_r", "").lower()
                folder.mkdir(exist_ok=True, parents=True)
                map_to_browse(
                    fits.open(path)[rec[0]].data,
                    latbins=latbins,
                    lonbins=lonbins,
                    outpath=folder / fn.lower(),
                    title=title + suffix,
                    cmap=cmap,
                    clip=clip,
                    normtype=norms.get(ptype, "standard"),
                    dpi_downscale_factor=dpi_downscale,
                    fontproperties=fontproperties,
                )


SETTINGS = {
    "fontsize": 4.5,
    "font": "TitilliumWeb-Regular.ttf",
    "slice_path": Path(__file__).parent.parent / "pipeline_data" / "deconv_slices",
    "map_path": Path(__file__).parent.parent / "bundle/data/maps",
    "browse_path": Path(__file__).parent.parent / "bundle/browse",
    "variable_meta": ("OLTSTMIN", "OLTSTMAX"),
    "max_abslat": 75,
    "ppd": 32,
    "dpi_downscale": 3,
    # NOTE: probably drop threads to 2 for 32 GB of RAM, and None for 16.
    # threads other than None, 2, or 4 is basically pointless because the
    # script parallelizes by channel and every channel (should) take the exact
    # same amount of iron.
    "threads": 4,
}

# browse product settings -- too annoying to set from CLI ever, just edit the
# file
CLIPS = {
    "STDEV": {"type": "centile", "value": (0, 99)},
    "DATMINUS": {"type": "centile", "value": (2, 98)},
    "LATSHIFT": {"type": "absolute", "value": (-15, 15)},
    "TEMP": {"type": "absolute", "value": (145, 290)},
    "TBMOD": {"type": "absolute", "value": (180, 305)},
}
CMAPS = {
    "LATSHIFT": ("Greys_r", "orange_teal"),
    "DATMINUS": ("inferno", "neon_r"),
    "TEMP": "moonbow",
    "STDEV": "vaportrail",
    "TBMOD": "sunset",
}
NORMS = {"LATSHIFT": "centered", "DATMINUS": "centered"}
FTYPES = ("base", "latshift", "browse", "datminus")
BROWSETYPES = ("TEMP", "STDEV", "LATSHIFT", "DATMINUS", "TBMOD")


def _assemble_channel_products(
    channel,
    browse_path,
    dpi_downscale,
    map_path,
    font,
    fontsize,
    slice_path,
    max_abslat,
    orbiter,
    ppd,
    tobrowse,
    tomake,
    variable_meta,
):
    mkwargs = {
        "orbiter": orbiter, "map_path": map_path, "channel": channel
    }
    watch = Stopwatch()
    watch.start()
    if "base" in tomake:
        print(f"****creating temp/stdev FITS for channel {channel}****")
        write_base_hdul(
            **mkwargs,
            slice_path=slice_path,
            max_abslat=max_abslat,
            variable_meta=variable_meta,
            ppd=ppd
        )
        print(f"done with temp/stdev ({watch.clickpeek()})")
    if "latshift" in tomake:
        print(f"****creating latshift FITS for channel {channel}****")
        write_latshift_hdul(**mkwargs, ppd=ppd)
        print(f"done with latshift ({watch.clickpeek()})")
    if "datminus" in tomake:
        print(f"****creating datminus FITS for channel {channel}****")
        write_datminus_hdul(**mkwargs, ppd=ppd)
        print(f"done with datminus ({watch.clickpeek()})")
    if "browse" in tomake:
        print(f"****creating browse products for channel {channel}****")
        make_browse_products(
            channel=channel,
            map_path=map_path,
            browse_path=browse_path,
            ptypes=tobrowse,
            fontproperties=mplf.FontProperties(
                fname=FONT_PATH / font, size=fontsize
            ),
            cmaps=CMAPS,
            norms=NORMS,
            clips=CLIPS,
            dpi_downscale=dpi_downscale,
            orbiter=orbiter
        )
        f"done with browse products ({watch.clickpeek()})"
    print(f"done with channel {channel} ({round(watch.total, 2)}s)")


def _run(
    orbiter,
    channel,
    make,
    browse_types,
    *,
    slice_path,
    map_path,
    browse_path,
    variable_meta,
    max_abslat,
    ppd,
    fontsize,
    font,
    dpi_downscale,
    threads
):
    tomake = FTYPES if make == "all" else tuple(map(str.lower, listify(make)))
    if not set(FTYPES).issuperset(tomake):
        raise ValueError(f"bad filetype specification {make}")
    # TODO, maybe: also check this for browse extensions but :shrug:
    if browse_types == "all":
        tobrowse = BROWSETYPES
    else:
        tobrowse = tuple(map(str.upper, listify(browse_types)))
    channels = (1, 2, 3, 4) if channel == "all" else listify(channel)
    pool = MaybePool(threads)
    kwargs = {
        "browse_path": browse_path,
        "dpi_downscale": dpi_downscale,
        "map_path": map_path,
        "font": font,
        "fontsize": fontsize,
        "slice_path": slice_path,
        "max_abslat": max_abslat,
        "orbiter": orbiter,
        "ppd": ppd,
        "tobrowse": tobrowse,
        "tomake": tomake,
        "variable_meta": variable_meta
    }
    bigwatch = Stopwatch()
    bigwatch.start()
    pool.map(
        _assemble_channel_products,
        [
            {'key': c, 'kwargs': {'channel': c} | kwargs}
            for c in channels
        ]
    )
    pool.close()
    pool.join()
    for k, v in pool.get().items():
        if isinstance(v, Exception):
            print(f"**** failed on channel {k} ****")
            rprint(exc_report(v))
            print("\n")
    print(f"****all done ({bigwatch.clickpeek()})****")


MapType = Literal["temp", "latshift", "datminus"]
Channel = Literal[1, 2, 3, 4]
MakeType = MapType | Literal["browse"]
BrowseType = Literal[MapType] | Literal["tbmod", "stdev"]


def run(
    *,
    orbiter: Literal["ce1", "ce2"],
    channel: Literal["all"] | tuple[Channel, ...] = "all",
    make: Literal["all"] | tuple[MakeType, ...] = "all",
    browse_types: Literal["all"] | tuple[BrowseType, ...] = "all",
    **kwargs
):
    _run(
        orbiter, channel, make, browse_types, **(SETTINGS | kwargs)
    )


if __name__ == "__main__":
    fire.Fire(run)
    # run(orbiter="ce2", channel=(2,), make=("latshift",))


