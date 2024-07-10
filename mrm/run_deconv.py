from itertools import product
from pathlib import Path
import pickle
import re
from typing import Collection, Literal, Union

import pandas as pd
from astropy.io import fits
from dustgoggles.dynamic import exc_report
from dustgoggles.structures import listify
import fire
from hostess.monitors import Stopwatch
import numpy as np
from rich import print as rprint

from mrm.processing.deconvolution import DeconvManager


def mgr_to_fn(mgr: DeconvManager):
    ltst = map(lambda r: r * 24, mgr.ltst_bounds)
    pnames = {
        n: "_".join(map(lambda coord: str(round(coord, 1)), a))
        for a, n in zip(
            (mgr.latrange, mgr.lonrange, ltst), ("lat", "lon", "ltst")
        )
    }
    fn = (
        (
            f"{mgr.orbiter}_t{mgr.channel}_ltst_{pnames['ltst']}_"
            f"lat_{pnames['lat']}_lon_{pnames['lon']}_deconv"
        )
        .replace(".0", "")
        .replace("-", "m")
    )
    return re.sub(r"[)(, '.]", "_", fn)


def write_deconv_fits(mgr, path):
    # noinspection PyTypeChecker
    hdus = [
        fits.PrimaryHDU(),
        fits.ImageHDU(np.flip(mgr.temp, axis=0), name=f"T{mgr.channel}"),
        fits.ImageHDU(np.flip(mgr.std, axis=0), name=f"T{mgr.channel}STDEV"),
        fits.ImageHDU(np.flip(mgr.fullbins["lat"]), name="LATBINS"),
        fits.ImageHDU(mgr.fullbins["lon"], name="LONBINS"),
    ]
    hdul = fits.HDUList(hdus)
    hdul["PRIMARY"].header["OLTSTMIN"] = mgr.ltst_bounds[0]
    hdul["PRIMARY"].header["OLTSTMAX"] = mgr.ltst_bounds[1]
    hdul["PRIMARY"].header["CHANNEL"] = mgr.channel
    hdul["PRIMARY"].header["DMINGAIN"] = mgr.mapped_gain_cutoff
    hdul["PRIMARY"].header["SIGMA"] = mgr.sigma
    hdul["PRIMARY"].header["ALTAZRES"] = mgr.altaz_resolution
    hdul["PRIMARY"].header["MINGAIN"] = mgr.gain_cutoff
    hdul["PRIMARY"].header["MAINBEAM"] = mgr.mainbeam_only
    hdul["PRIMARY"].header["PPD"] = round(1 / mgr.latlon_resolution)
    hdul["PRIMARY"].header["ORBITER"] = mgr.orbiter
    if mgr.extended_info is True:
        hdul.append(
            fits.ImageHDU(np.flip(mgr.scounts, axis=0), name="SCOUNTS")
        )
        orbs = map(
            lambda l: np.array(sorted(l)).astype(np.uint16), mgr.orbs.ravel()
        )
        import pyarrow as pa
        from pyarrow import parquet as pq
        pq.write_table(
            pa.Table.from_pandas(
                pd.DataFrame({'orbits': list(orbs)}), preserve_index=False
            ),
            path.parent / f"{path.stem}_orbtab.parquet"
        )
    hdul.writeto(path, overwrite=True)


def make_deconv_slice(
    outpath, n_threads, skip_existing, dump_failures, **manager_kwargs
):
    manager = DeconvManager(**manager_kwargs)
    path = outpath / f"{mgr_to_fn(manager)}.fits"
    if skip_existing is True and path.exists():
        print(f"{outpath / mgr_to_fn(manager)} exists, skipping...")
        return
    elif needs_lockfile := not path.exists():
        path.touch()
    print(manager)
    try:
        manager.prep()
        print(len(manager.table))
        manager.run(n_threads)
        if manager._done is True:
            write_deconv_fits(manager, path)
        if dump_failures is True:
            with open(f"{mgr_to_fn(manager)}.pkl", "wb") as stream:
                pickle.dump(manager.failures, stream)
    finally:
        if manager._done is False and needs_lockfile is True:
            path.unlink(missing_ok=True)


def run(
    orbiter: Literal["ce1", "ce2", "all"] = "all",
    channel: Literal[1, 2, 3, 4, "all"] = "all",
    tbin_size: int = 2,
    tbin_ix: Union[Literal["all"], int, Collection[int]] = "all",
    ppd=32,
    skip_existing=False,
    **kwargs
):
    kwargs = SETTINGS | kwargs
    kwargs['outpath'] = Path(kwargs['outpath'])
    kwargs['outpath'].mkdir(exist_ok=True)
    kwargs |= {"latlon_resolution": 1/ppd}
    if 24 % tbin_size != 0:
        raise ValueError("tbin_size must be a divisor of 24.")
    tbin_ix = listify(tbin_ix)
    watch = Stopwatch(digits=2)
    watch.start()
    channels = (channel,) if channel != "all" else (1, 2, 3, 4)
    orbiters = (orbiter,) if orbiter != "all" else ("ce2", "ce1")
    tbins = tuple(range(0, int(24 / tbin_size)))
    for orb, chan, (tix, tbin) in product(
        orbiters, channels, enumerate(tbins)
    ):
        if tbin_ix != ["all"] and tix not in tbin_ix:
            continue
        kwargs['orbiter'], kwargs['channel'] = orb, chan
        lbounds = (
            tbin * tbin_size / 24, (tbin * tbin_size + tbin_size) / 24
        )
        try:
            make_deconv_slice(
                **kwargs, ltst_bounds=lbounds, skip_existing=skip_existing
            )
        except KeyboardInterrupt:
            raise
        except Exception as ex:
            rprint("[bold red]!!!! deconvolution failed !!!![/bold red]\n")
            rprint(exc_report(ex))
            print("\n")
        print(f"\n****{watch.clickpeek()}****")


SETTINGS = {
    "altaz_resolution": 0.05,
    "sigma": 3,
    "latrange": (-75, 75),
    "lonrange": (-180, 180),
    "unflagged": True,
    "mapped_gain_cutoff": -9,
    "gain_cutoff": -24,
    "max_pattern_angle": 70,  # should not matter unless really small
    "mainbeam_only": True,
    "tilesize": 4,
    "n_threads": 4,
    "base_tablepath": (
        Path(__file__).parent.parent
        / "pipeline_data/concatenated_table/ce_mrm_concatenated.parquet"
    ),
    "outpath": (
        Path(__file__).parent.parent / "pipeline_data/deconv_slices"
    ),
    "extended_info": False,
    "dump_failures": False,
}


if __name__ == "__main__":
    fire.Fire(run)