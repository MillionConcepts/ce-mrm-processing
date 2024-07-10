import atexit
import time
import warnings
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from functools import partial, wraps
from inspect import getfullargspec
from itertools import chain, product
from numbers import Number, Real
from pathlib import Path
from typing import Callable, Collection, Literal, Optional

import pandas as pd
from cytoolz import curry, get, merge, valmap, first
from dustgoggles.structures import HashDict, MaybePool
import fast_histogram as fh
from hostess.monitors import Stopwatch
from marslab.imgops.imgutils import ravel_valid
from moonbow import div0  # we love dividing by 0
from more_itertools import windowed
import numba as nb
import numpy as np
import pyarrow as pa
from pyarrow import compute as pc, parquet as pq
import scipy.ndimage as ndi

from mrm.processing.numba_cheats import njit_dyncached
from mrm.processing.projection import qasin, qatan2, qcos, qsin
from mrm.shared.console import print_inline
from mrm.shared.constants import ANTENNA_PARAMS, LUNAR_RADIUS_KM, RADIANS
from mrm.shared.mrmtypes import MRMChannel
from scipy.stats import binned_statistic_2d


def db2power(db: np.ndarray | Number) -> np.ndarray | Number:
    """Convert dB to normalized power."""
    return 10 ** (db / 10)


def power2db(power: np.ndarray | Number) -> np.ndarray | Number:
    """Convert normalized power to dB."""
    return 10 * np.log10(power)


def biased_std(
    weighted_values: np.ndarray,
    weighted_squares: np.ndarray,
    weights: np.ndarray
) -> np.ndarray:
    """
    Calculate biased estimate of weighted standard deviation.

    Args:
        weighted_values: sum of weighted values
        weighted_squares: sum of weighted squares of values
        weights: sum of weights

    Returns:
        Biased estimate of weighted standard deviation.
    """
    return np.sqrt(
        (weighted_squares * weights - weighted_values**2) / weights**2
    )


def _crop_to_bound(
    bound_ix: int, gain: np.ndarray, theta: np.ndarray
) -> tuple[Real, np.ndarray, np.ndarray]:
    """
    Perform a square crop on `gain` and `theta`.

    Args:
        bound_ix: offset from 'center' of arrays at which to crop them.
        gain: a 2D array to crop; in the antenna pattern pipeline, it
            represents gain.
        theta: a 2D array to crop and retrieve a bounding value from.
            In the antenna pattern pipeline, it represents magnitude of alt/az
            offset from boresight. The valid range of the model is a circle,
            so cropping to its rectangular extent leaves invalid values in the
            corners; this function also grabs a bounding value from `theta` to
            later set them to NaN. Note that this assumes that `theta` is
            radially symmetrical.

    Returns:
         bound_theta: `theta` at (center index, center index + `bound_ix`)
         cropped_gain: cropped copy of `gain`
         cropped_theta: cropped copy of `theta`
    """
    center_ix = int(len(theta) / 2)
    bound_theta = theta[center_ix, center_ix + bound_ix]
    # crop arrays to valid range
    gain = gain[
        center_ix - bound_ix : center_ix + bound_ix + 1,
        center_ix - bound_ix : center_ix + bound_ix + 1,
    ]
    theta = theta[
        center_ix - bound_ix : center_ix + bound_ix + 1,
        center_ix - bound_ix : center_ix + bound_ix + 1,
    ]
    return bound_theta, gain, theta


def ce2_antenna_pattern(
    theta: nb.float32[:, :],
    mainbeamwidth: nb.float32,
    sigma: nb.float32,
    m_power: nb.float32,
    sidepeak: nb.float32,
    s_by_m: nb.float32,
    mainbeam_only: bool = False,
    gain_cutoff: Optional[float] = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """

    Create power, gain, and angular offset arrays defining an MRM antenna
    pattern. All arguments other than `theta`, `mainbeam_only`, and
    `gain_cutoff` are model parameters; their canonical values per MRM channel
    are defined in `ANTENNA_PARAMS`.

    Args:
        theta: Array containing magnitudes of alt/az offsets from boresight.
            Units are degrees.
        mainbeamwidth: Extent of modeled main beam. Units are degrees.
        sigma: Spread term. Used to construct parameters of Gaussian
            distributions.
        m_power: Mainbeam power term. Unitless.
        sidepeak: Sidelobe power term. Unitless.
        s_by_m: Sidelobe spread term. Unitless.
        mainbeam_only: If True, only model mainbeam.
        gain_cutoff: If not None, crop returned arrays to the largest square
            such that all gain values outside the square are less than
            `gain_cutoff`.

    Returns:
        power: ndarray of pattern weights expressed as normalized power
        gain: ndarray of pattern weights expressed as gain / dB, referenced
            to maximum instantaneous power
        pattern_theta: ndarray giving alt/az offset from boresight in degrees;
            this is the input argument `theta` cropped to the valid extent of
            the model.

    Notes:
        Although the first MRM sidelobe is highly asymmetrical, this function
        models it as radially symmetrical, basically giving its average value
        at all rotational positions. The data that would be required to
        determine beam orientation are not available, so it is impossible, in
        practice, to do better than this. Due to this uncertainty, adding the
        sidelobe term makes the accuracy of the deconvolution slightly worse
        at the cost of many processor cycles. As such, the default
        deconvolution pipeline configuration sets `mainbeam_only` to `True`.
        Approximating the first sidelobe may still be useful for some purposes,
        however, so we have retained this portion of the model as an option.

        For details on assumptions about antenna patterns, see:

        Wang et al. (2010). Calibration and brightness temperature algorithm
        of CE-1 Lunar Microwave Sounder (CELMS). Science China Earth Sciences,
        53, 392–1406. https://doi.org/10.1007/s11430-010-4008-x
    """
    halfmb = mainbeamwidth / 2
    m_diff = m_power - sidepeak
    center_ix = int(len(theta) / 2)  # nominal center index of array
    mainbeam_mask = theta < halfmb  # indices at which main beam is defined
    # 'canvas' for constructing model
    gain = np.zeros_like(theta)
    if mainbeam_only is False:
        # construct first sidelobe
        sidelobe_angle = (theta - halfmb) / s_by_m
        # subtractive model: inner_annulus and outer_annulus
        # express how the first sidelobe falls from peak gain
        # (`sidepeak`)
        inner_annulus = m_diff * np.exp(
            -((sidelobe_angle / sigma) ** nb.f4(2))
        )
        outer_annulus = m_diff * np.exp(
            -(((sidelobe_angle - mainbeamwidth) / sigma) ** nb.f4(2))
        )
        gain += outer_annulus + inner_annulus + sidepeak
    # construct main beam
    center_theta = theta[mainbeam_mask]
    # this term doesn't do a ton but it's in the model
    boostpeak = m_power * np.exp(-(((center_theta + halfmb) / sigma) ** 2))
    mainbeam = m_power * np.exp(-(((center_theta - halfmb) / sigma) ** 2))
    gain[mainbeam_mask] = boostpeak + mainbeam
    # crop to valid region
    gain_ratio_cutoff = 0.01  # arbitrary
    if mainbeam_only is False:
        # The subtractive annulus functions hits limits and stop subtracting,
        # meaning that after the outer annulus hits its zero, the function
        # begins increasing again, and the 'sidepeak' value extends to
        # infinity as spurious background. we want to cut off after this
        # limit. function is radially symmetrical, so we just check along
        # the x axis at y center.
        # noinspection PyUnboundLocalVariable
        annulus_line = outer_annulus[center_ix, center_ix:]
        bound_ix = np.nonzero(
            np.abs(annulus_line - annulus_line.min()) < gain_ratio_cutoff
        )[0].min()
    else:
        # if we only care about the main beam, it's just its nominal bound:
        bound_ix = np.argmax(~mainbeam_mask[center_ix, center_ix:]) - 1
    bound_theta, gain, theta = _crop_to_bound(bound_ix, gain, theta)
    # optionally add gain as a boundary condition
    if gain_cutoff is not None:
        center_ix = int(len(theta) / 2)
        low = np.nonzero(gain[center_ix, center_ix:] < gain_cutoff)[0]
        if len(low) > 0:
            bound_ix = (
                np.split(low, np.nonzero(np.diff(low) != 1)[0] + 1)[-1].min()
                - 1
            )
            bound_theta, gain, theta = _crop_to_bound(bound_ix, gain, theta)
    # set corners to -inf dB
    gain[theta > bound_theta] = -np.inf
    # convert from dB to normalized power
    return db2power(gain - gain.max()), gain, theta


def make_pattern(
    channel: MRMChannel,
    max_angle_altaz: float = 28,
    pattern_resolution_altaz: float = 0.01,
    mainbeam_only: bool = False,
    gain_cutoff: Optional[float] = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Handler function for antenna pattern modeling.

    Args:
        channel: which channel's pattern to model
        max_angle_altaz: alt/az cutoff for pattern, in degrees. Should not
        pattern_resolution_altaz: resolution of pattern, in degrees.
        mainbeam_only: model only the mainbeam?
        gain_cutoff: if not None, square-crop the pattern at the point that
            all values fall under `gain_cutoff`.

    Returns:
        power: ndarray of pattern weights expressed as normalized power
        gain: ndarray of pattern weights expressed as gain / dB, referenced
            to maximum instantaneous power
        theta: ndarray giving alt/az offset from boresight in degrees
    """
    # define beam pattern in altitude, azimuth.
    # azimuth/elevation axis is in degrees.
    azax = np.arange(
        -max_angle_altaz,
        max_angle_altaz + pattern_resolution_altaz,
        pattern_resolution_altaz,
        dtype="f4",
    )
    # deflection from boresight in degrees
    theta = sum(map(lambda ax: ax**2, np.meshgrid(azax, azax))) ** 0.5
    return ce2_antenna_pattern(
        theta,
        **valmap(np.float32, ANTENNA_PARAMS[channel]),
        gain_cutoff=gain_cutoff,
        mainbeam_only=mainbeam_only,
    )


# noinspection PyUnresolvedReferences
def azel_to_latlon_source(
    orbalt: nb.float32,
    lat_0: nb.float32,
    lon_0: nb.float32,
    y_over_rad: nb.float32[:],
    x_over_rad: nb.float32[:],
    scaled_rad: nb.float32[:]
) -> tuple[nb.float32[:], nb.float32[:]]:
    """
    Source code for accelerated version of reverse orthographic projection
    function. Projects instrument footprint, measured in az/el relative to
    boresight, to lat/lon on planetary surface. The compiled version provides
    an order-of-magnitude speedup over the unaccelerated version (more at
    lower latitudes).

    Approximates the surface as locally flat and assumes the instrument is
    looking perfectly nadir. This allows us to treat the offset from center
    as a right triangle, so linear surface extent in rectangular coordinates
    is simply tan(angular offset from boresight / 2) * orbital altitude * 2.

    Then we assume that the target body is perfectly spherical and transform
    these rectangular coordinates into lat/lon using standard spherical trig
    formulae.

    Args:
        orbalt: orbital altitude [km]
        lat_0: latitude of boresight intercept point [deg]
        lon_0: longitude of boresight intercept point [deg]
        y_over_rad: ratios of y-axis vector components to vector magnitudes
            [unitless]
        x_over_rad: ratios of x-axis vector components to vector magnitudes
            [unitless]
        scaled_rad: vector magnitudes [km]

    Returns:
        lat: array of latitude values [deg]
        lon: array of longitude values [deg]

    Notes:
        1. In the deconvolution pipeline, this serves as a prototype for
            compiled functions partially evaluated with precomputed values of
            `y_over_rad`, `x_over_rad`, and `scaled_rad`; this function itself
            is never called.
        2. For performance reasons, this function works with and returns
            "raveled" (1-D) arrays and is unaware of the shape of the
             original 2-D arrays. Transforming them into 2-D arrays, if
            required, is the responsibility of the caller.
    """
    lat_0 = lat_0 * RADIANS
    lon_0 = lon_0 * RADIANS
    c = scaled_rad * orbalt
    cos_c, sin_c = qcos(c), qsin(c)
    # arcsin approximation by identity function is not good here,
    # specifically, it breaks at higher lat
    lat = qasin(
        cos_c * np.sin(lat_0)
        + y_over_rad * sin_c * np.cos(lat_0)
    )
    lon = lon_0 + qatan2(
        x_over_rad * sin_c,
        cos_c * np.cos(lat_0) - y_over_rad * sin_c * np.cos(lat_0)
    )
    return lat / RADIANS, lon / RADIANS


AZEL_TO_LATLON_SIG = nb.types.UniTuple(nb.float32[:], 2)(
    nb.float32,
    nb.float32,
    nb.float32,
    nb.float32[:],
    nb.float32[:],
    nb.float32[:]
)
"""numba Signature for compiled versions of azel_to_latlon"""


PIPELINE_CONSTANTS = (
    "ALTAZ_EXTENT_RAD",
    "EXTENT_CONSTANT",
    "NPIX",
    "POWER",
    "GAIN",
    "THETA",
    "SCALED_RAD",
    "X_OVER_RAD",
    "Y_OVER_RAD",
    "WEIGHTS",
)
"""
Names of precomputed constants for projection / deconvolution pipeline.

TODO: not all of these are actually required; remove those that are not.
"""

A2L_CONSTANTS = ("Y_OVER_RAD", "X_OVER_RAD", "SCALED_RAD")
"""
Names of precomputed constants intended to be bound in the enclosing scopes of 
compiled versions of azel_to_latlon.
"""


def prep_pipeline(
    channel: MRMChannel,
    pattern_resolution_altaz: float,
    max_pattern_angle: float = 80,
    mainbeam_only: bool = True,
    gain_cutoff: Optional[float] = None,
) -> dict[str, np.ndarray | Real]:
    """
    Compute antenna patterns and reverse orthographic projection constants for
    deconvolution pipeline. Intended primarily be called by  `DeconvManager`,
    but also usable on its own.

    Args:
        channel: MRM channel
        pattern_resolution_altaz: alt/az resolution for antenna pattern [deg]
        max_pattern_angle: crop pattern at this alt/az magnitude [deg];
            should not matter unless very small
        mainbeam_only: if True, model only the mainbeam; otherwise, also model
            the first sidelobe.
        gain_cutoff: if not None, crop the model to the square extent at which
            no point outside that extent has dB >= `gain_cutoff`.

    Returns:
        dict of pipeline constants.
    """
    # make antenna pattern and angular offset backplane
    power, gain, theta = make_pattern(
        max_angle_altaz=max_pattern_angle,
        channel=channel,
        pattern_resolution_altaz=pattern_resolution_altaz,
        mainbeam_only=mainbeam_only,
        gain_cutoff=gain_cutoff,
    )
    # make constants for azel_to_latlon
    iz = np.ones(power.shape, dtype=np.float32)
    iz[power == 0] = np.nan
    npix = np.float32(power.shape[0])
    altaz_extent_rad = np.float32(
        np.radians(theta[int(theta.shape[0] / 2), 0])
    )
    ix_y, ix_x = np.indices(power.shape)
    ix_y, ix_x = map(lambda c: c - npix / 2 + 0.5, (ix_y, ix_x))
    ix_y, ix_x = map(lambda c: c * iz, (ix_y, ix_x))
    rix_y, rix_x = map(
        lambda c: ravel_valid(c.astype(np.float32)), (ix_y, ix_x)
    )
    weights = ravel_valid(power * iz)
    # precomputed terms related to distance from boresight
    rix_rad = (rix_x**2 + rix_y**2) ** 0.5 + 1e-9  # prevent annoying NaN
    extent_constant = np.tan(altaz_extent_rad) * np.float32(2) / npix
    scaled_rad = extent_constant * rix_rad / LUNAR_RADIUS_KM
    y_over_rad, x_over_rad = rix_y / rix_rad, rix_x / rix_rad
    return {
        k.upper(): v
        for k, v in locals().items()
        if k.upper() in PIPELINE_CONSTANTS
    }


def pq_between(
    field: str, bounds: tuple[Real, Real]
) -> list[tuple[str, Literal[">=", "<"], Real]]:
    """
    Shorthand function that makes an object that can be passed to
    `pq.read_table()` as its `filters` argument to specify
    `bounds[0] <= value of `field` < bounds[1]`.

    Args:
        field: name of field of table schema on which to filter
        bounds: tuple giving bounds of filter interval

    Returns:
        Suitable `filters` for `pq.read_table()`.
    """
    return [(field, ">=", bounds[0]), (field, "<", bounds[1])]


def get_pattern_spill(
    max_orbalt: Real,
    max_abslat: Real,
    a2l: Callable,
    a2l_constants: dict[str, np.ndarray | Real]
) -> dict[str, np.float32]:
    """
    For a particular antenna pattern definition (defined via precomputed
    projection constants), orbital altitude, and absolute latitude, what is
    the maximum lat/lon offset between the boresight and the edge of the
    antenna pattern? Used to add worst-case "spill" bins to tiles initially
    defined by point-sample / boresight coordinates.
    """
    # (note that this is not longitude-variant)
    testlat, testlon = a2l(
        np.float32(max_orbalt),
        np.float32(max_abslat),
        np.float32(30),  # arbitrary, literally does not matter
        *a2l_constants
    )
    return {"lat": np.ptp(testlat) / 2, "lon": np.ptp(testlon) / 2}


def make_full_grid(
    grid_resolution: Real,
    lat_bounds: tuple[Real, Real],
    lon_bounds: tuple[Real, Real],
    spill: dict[Literal["lat", "lon"], Real],
    pin_bin_edges: bool = True
) -> dict[str, np.ndarray]:
    """
    Produce lat/lon bins (equivalently, discrete axes) for an equirectangular
    projection grid, adding "spill bins" to account for extent of effective
    measurement area from nominal geolocated values of sample population.

    Args:
        grid_resolution: side length of each square grid cell [deg]
        lat_bounds: base min/max grid latitude [deg]
        lon_bounds: base min/max grid longitude [deg, -180 to 180 system]
        spill: dict giving maximum possible 'spill' along each of latitude and
            longitude axes in degrees.
        pin_bin_edges: if True, shift bins so that the leftmost lat/lon bins
            are at integer offsets from lat/lon

    Returns:
        dict like: {'lat': lat_bins, 'lon': lon_bins}, where lat_bins and
            lon_bins are 1-D ndarrays.

    Notes:
        1. This function, by intent, will happily produce out-of-range
        longitude bins. This is an intentional optimization for the pipeline;
        it allows it to skip bounds checks and correctly "wrap" bins at the
        end.
        2. Returned array values should represent left edges in all cases.
    """
    bins = {
        "lat": np.arange(
            # TODO: this shouldn't need to be symmetrical in lat; there's some
            #  array initialization bug that makes it throw errors when not.
            -max(lat_bounds) - spill["lat"],
            max(lat_bounds) + spill["lat"],
            grid_resolution,
            dtype="f4",
        ),
        "lon": np.arange(
            lon_bounds[0] - spill["lon"],
            lon_bounds[1] + spill["lon"],
            grid_resolution,
            dtype="f4",
        ),
    }
    if pin_bin_edges is True:
        bins["lon"] -= bins["lon"].min() % 1
        bins["lat"] -= bins["lat"].min() % 1
    return bins


# define a smaller grid with aligned bins
def get_tile_bins(full_bins, bounds, max_altitude, a2l, a2l_constants):
    tspill = get_pattern_spill(
        max_altitude,
        np.abs(np.array(bounds["lat"], dtype=np.float32)).max(),
        a2l,
        a2l_constants
    )
    binslice = {}
    for ax in ("lat", "lon"):
        binrange = []
        bottom = np.nonzero(full_bins[ax] < bounds[ax][0] - tspill[ax])[0]
        if (len(bottom) == 0) or (bottom[-1] == 0):
            binrange.append(0)
        else:
            binrange.append(bottom[-1] - 1)
        top = np.nonzero(full_bins[ax] > bounds[ax][1] + tspill[ax])[0]
        if (len(top) == 0) or (top[0] == len(full_bins) - 2):
            binrange.append(len(full_bins[ax]) - 1)
        else:
            binrange.append(top[0] + 1)
        binslice[ax] = tuple(binrange)
    bdicts = [
        {ax: full_bins[ax][slice(*binslice[ax])], f"{ax}_bins": binslice[ax]}
        for ax in ("lat", "lon")
    ]
    return merge(bdicts)


# NOTE: if we ever make this extensible to irregular grids,
#  we will have to switch back from fast-histogram to sumbins().
# noinspection PyUnboundLocalVariable
def integrate_tile(
    tframe,
    tbins,
    a2l,
    sigma,
    weights,
    a2l_constants,
    extended: bool = False,
):
    num = np.zeros((len(tbins["lat"]), len(tbins["lon"])), dtype="f4")
    num2 = num.copy()
    denom = num.copy()
    if extended is True:
        scounts = np.zeros(num.shape, dtype="uint32")
        # orbs = np.full(num.shape, None, dtype="object")
        # for i, j in np.ndindex(orbs.shape):
        #     orbs[i, j] = set()
    latrange = tbins["lat"].min(), tbins["lat"].max()
    n_latbins = tbins["lat"].size
    lonrange = tbins["lon"].min(), tbins["lon"].max()
    n_lonbins = tbins["lon"].size
    for i, row in tframe.iterrows():
        lat, lon = a2l(row["d"], row["lat"], row["lon"], *a2l_constants)
        integrated = fh.histogram2d(
            lat, lon, (n_latbins, n_lonbins), (latrange, lonrange), weights
        )
        if sigma is not None:
            integrated = ndi.gaussian_filter(integrated, sigma=sigma)
        num += integrated * row["t"]
        num2 += integrated * row["t"] ** 2
        denom += integrated
        if extended is True:
            scounts[integrated >= 7] += 1
            # orb = int(row['orbit'])
            # for val in orbs[integrated > 0].ravel():
            #     val.add(orb)
    output = {
        "weighted_values": num, "weighted_squares": num2, "counts": denom
    }
    if extended is False:
        return output
    return output | {'scounts': scounts}


def load_slice(
    tablepath: Path,
    geo_bounds: dict[str, tuple[float, float]],
    channel: int,
    orbiter: str = None,
    ltst_bounds=None,
    unflagged=True,
    exclude_orbits=None,
    presliced=False,
) -> pa.Table:
    """
    Load a 'slice' of one of the concatenated data parquet files (or a 'fast'
    table generated from one). Used in the tile-processing workflow.
    """
    filters = list(
        chain(*[pq_between(ax, geo_bounds[ax]) for ax in ("lat", "lon")])
    )
    if presliced is False:
        filters.append(("orbiter", "=", orbiter))
        if unflagged is True:
            filters.append(("flag", "=", 0))
        if ltst_bounds is None:
            tablespecs = [filters]
        elif min(ltst_bounds) >= 0:
            tablespecs = [filters + pq_between("ltst", ltst_bounds)]
        else:
            tablespecs = [
                filters + pq_between("ltst", p)
                for p in [(1 + min(ltst_bounds), 1), (0, max(ltst_bounds))]
            ]
    else:
        tablespecs = [filters]
    columns = ["d", "lat", "lon", "ltst", f"t{channel}", "orbit"]
    slicetables = [
        pq.read_table(tablepath, columns=columns, filters=s)
        for s in tablespecs
    ]
    concat = pa.concat_tables(slicetables)
    if exclude_orbits is not None:
        concat = concat.filter(
            pc.invert(pc.is_in(concat["orbit"], pa.array(exclude_orbits)))
        )
    return concat


def load_tileframe(geo_bounds, **slice_kwargs):
    tileframe = load_slice(geo_bounds=geo_bounds, **slice_kwargs).to_pandas()
    tileframe.columns = ["d", "lat", "lon", "ltst", "t", "orbit"]
    return tileframe


def process_tile(
    weights: np.ndarray,
    nb_cachename: str,
    geo_bounds,
    a2l_constants,
    bins: dict[str, np.ndarray],
    sigma: Optional[int] = 3,
    extended: bool = False,
    **slice_kwargs,
):
    prepstart = time.time()
    tframe = load_tileframe(geo_bounds, **slice_kwargs)
    if len(tframe) == 0:
        raise ValueError("This tile has no matching data.")
    a2l, _ = njit_dyncached(
        azel_to_latlon_source, nb_cachename, sig=AZEL_TO_LATLON_SIG
    )
    tbins = get_tile_bins(
        bins, geo_bounds, tframe["d"].max(), a2l, a2l_constants
    )
    preptime = time.time() - prepstart
    rec = {"tbins": tbins, "preptime": preptime} | slice_kwargs
    slicestart = time.time()
    integration = integrate_tile(
        tframe,
        tbins,
        a2l,
        sigma,
        weights,
        a2l_constants,
        extended
    )
    return rec | integration | {"slicetime": time.time() - slicestart}


def prep_fast_table(
    base_tablepath: Path,
    geo_bounds,
    channel,
    orbiter,
    ltst_bounds=None,
    unflagged=True,
    exclude_orbits=None,
):
    """
    make a parquet file optimized for fast loading in the deconv pipeline.
    """
    chopped = load_slice(
        base_tablepath,
        geo_bounds,
        channel,
        orbiter,
        ltst_bounds,
        unflagged,
        exclude_orbits,
    ).sort_by("lat")
    outpath = Path(
        f".temp/ftab{str(np.random.randint(10000)).zfill(4)}.parquet"
    )
    outpath.parent.mkdir(exist_ok=True, parents=True)
    pq.write_table(chopped, outpath, row_group_size=10000)
    return outpath


def cut_valid(arr):
    valid_y, valid_x = np.nonzero(np.isfinite(arr))
    vslice = (
        slice(valid_y.min(), valid_y.max()),
        slice(valid_x.min(), valid_x.max()),
    )
    return {"y": valid_y, "x": valid_x}, vslice, arr[*vslice]


def maybedelete(path):
    if path.exists():
        path.unlink()


@curry
def only_when(method, attr, msg=None):
    if msg is None:
        msg = f"{method.__name__}() only valid when {attr.strip('_')}."

    @wraps(method)
    def run_only_when(obj, *args, **kwargs):
        if getattr(obj, attr) is False:
            raise ValueError(msg)
        return method(obj, *args, **kwargs)

    return run_only_when


class DeconvManager:

    def __init__(
        self,
        *,
        channel: Literal[1, 2, 3, 4],
        latrange: tuple[float, float],
        lonrange: tuple[float, float],
        base_tablepath: Path,
        ltst_bounds: Optional[tuple[float, float]] = None,
        exclude_orbits: Optional[Collection[int]] = None,
        unflagged: bool = True,
        altaz_resolution: float = 0.1,
        latlon_resolution: float = 1 / 64,
        sigma: Optional[int] = None,
        gain_cutoff: Optional[float] = None,
        max_pattern_angle: float = 70,
        mainbeam_only: bool = False,
        clean_tempfiles: bool = True,
        orbiter: Literal["ce1", "ce2"] = "ce2",
        tilesize: int = 4,
        pin_bin_edges: bool = True,
        mapped_gain_cutoff: Optional[float] = -9,
        extended_info: bool = False
    ):
        with self._setmode():
            self.channel, self.latrange = channel, latrange
            self.lonrange, self.altaz_resolution = lonrange, altaz_resolution
            self.latlon_resolution, self.sigma = latlon_resolution, sigma
            self.pin_bin_edges = pin_bin_edges
            self.gain_cutoff, self.mainbeam_only = gain_cutoff, mainbeam_only
            self.max_pattern_angle, self.orbiter = max_pattern_angle, orbiter
            self.clean_tempfiles, self.sigma = clean_tempfiles, sigma
            self.exclude_orbits, self.unflagged = exclude_orbits, unflagged
            self.tilesize, self.ltst_bounds = tilesize, ltst_bounds
            self.base_tablepath = base_tablepath
            self.mapped_gain_cutoff = mapped_gain_cutoff
            self.extended_info = extended_info

    def prep(self):
        with self._setmode():
            self.constants = prep_pipeline(
                pattern_resolution_altaz=self.altaz_resolution,
                max_pattern_angle=self.max_pattern_angle,
                channel=self.channel,
                mainbeam_only=self.mainbeam_only,
                gain_cutoff=self.gain_cutoff,
            )
            self.a2l, self.nb_cachename = njit_dyncached(
                azel_to_latlon_source,
                globals_={'RADIANS': RADIANS} | globals(),
                sig=AZEL_TO_LATLON_SIG,
            )
            # TODO, maybe: we should not really be restricting to samples
            #  that fall inside of the box. Most correctly, we should add
            #  spill for sample selection as well. However, these do not create
            #  issues for global maps.
            self.fast_tablepath = prep_fast_table(
                self.base_tablepath,
                {"lat": self.latrange, "lon": self.lonrange},
                self.channel,
                self.orbiter,
                self.ltst_bounds,
                self.unflagged,
                exclude_orbits=self.exclude_orbits,
            )
            if self.clean_tempfiles is True:
                atexit.register(partial(maybedelete, self.fast_tablepath))
            self.table = pq.read_table(self.fast_tablepath)
            if len(self.table) == 0:
                raise ValueError(
                    f"No data points in {self.base_tablepath.name} match "
                    f"requested parameters."
                )
            self.max_orbalt = np.float32(pc.max(self.table["d"]).as_py())
            self.spill = get_pattern_spill(
                self.max_orbalt,
                np.float32(pc.max(pc.abs(self.table["lat"])).as_py()),
                self.a2l,
                get(list(A2L_CONSTANTS), self.constants)
            )
            self.compbins = make_full_grid(
                self.latlon_resolution,
                self.latrange,
                self.lonrange,
                self.spill,
                self.pin_bin_edges,
            )
            self._make_argrecs()
            self._prepped = True

    @only_when(attr="_setting", msg="do not call this function directly")
    def _make_argrecs(self):
        self._make_tilespecs()
        basekwargs = {
            "presliced": True,
            "tablepath": self.fast_tablepath,
            "a2l_constants": get(list(A2L_CONSTANTS), self.constants),
            "weights": self.constants["WEIGHTS"],
            "bins": self.compbins,
            "channel": self.channel,
            "sigma": self.sigma,
            "nb_cachename": self.nb_cachename,
            "extended": self.extended_info,
        }
        self.argrecs = [
            {"kwargs": basekwargs | {"geo_bounds": tile}}
            for i, tile in enumerate(self.tilespecs)
        ]

    # TODO: something very weird is wrong with certain tilesize
    #  settings (like 3)
    def _geochunk(self, which):
        bounds = getattr(self, f"{which}range")
        return tuple(
            windowed(range(bounds[0], bounds[1] + 1, self.tilesize), 2)
        )

    @contextmanager
    def _setmode(self):
        self._setting = True
        try:
            yield []
        finally:
            self._setting = False

    @only_when(attr="_setting", msg="do not call this function directly")
    def _make_tilespecs(self):
        # noinspection PyTypeChecker
        self.chunks = self._geochunk("lon"), self._geochunk("lat")
        self.tilespecs = [
            {
                ax: (min(p[i]), max(p[i]))
                for ax, i in zip(("lon", "lat"), (0, 1))
            }
            for p in product(*self.chunks)
        ]

    def reset(self):
        with self._setmode():
            if self._running is True:
                self.pool.terminate()
            if self.clean_tempfiles is True and hasattr(self, "tpath"):
                maybedelete(self.fast_tablepath)
            for attr in self.prepattrs + self.runattrs:
                if hasattr(self, attr):
                    self.__delattr__(attr)
            for state in self.states:
                setattr(self, state, False)

    def __setattr__(self, attr, value):
        if (attr == "_setting") or self._setting is True:
            return super().__setattr__(attr, value)
        if attr in self.states + self.runattrs:
            raise ValueError("cannot set this attribute directly.")
        if (self._running or self._armed) is True:
            raise ValueError(
                "cannot set this attribute when armed, during run, "
                "or prior to run cleanup. Run reset() first."
            )
        self._prepped, self._readied = False, False
        self._setting = True
        for prepattr in self.prepattrs:
            if hasattr(self, prepattr):
                self.__delattr__(prepattr)
        super().__setattr__(attr, value)
        self._setting = False

    @only_when(attr="_running")
    def ready_results(self):
        return {k for k, v in self.pool.results_ready().items() if v is True}

    @only_when(attr="_done")
    def _make_backplanes(self):
        self.backplanes = {}
        self.backplanes["lat_full"] = np.tile(
            self.fullbins["lat"], (self.fullbins["lon"].size, 1)
        ).T
        self.backplanes["lon_full"] = np.tile(
            self.fullbins["lon"], (self.fullbins["lat"].size, 1)
        )
        self.backplanes["lat"] = self.backplanes["lat_full"][*self.vslice]
        self.backplanes["lon"] = self.backplanes["lon_full"][*self.vslice]

    @only_when(attr="_done")
    def _rectify_coordinates(self):
        outofbound_recs, okbins = [], set(range(self.compbins["lon"].size))
        for side, sign, op in zip(
            ("low", "high"), (-1, 1), (np.less, np.greater_equal)
        ):
            outbins = np.nonzero(op(self.compbins["lon"], 180 * sign))[0]
            if outbins.size == 0:
                continue
            okbins.difference_update(outbins)
            wraprec = {
                "lon": self.compbins["lon"][outbins] - 360 * sign,
                "outbins": outbins,
            }
            outofbound_recs.append(wraprec)
        if len(outofbound_recs) == 0:
            self.fullbins = self.compbins.copy()
            return
        splicelon = np.hstack(
            [*[w["lon"] for w in outofbound_recs], self.compbins["lon"]]
        )
        lonbins = np.arange(
            max(splicelon.min(), -180),
            min(splicelon.max(), 180 - self.latlon_resolution)
            + self.latlon_resolution,
            self.latlon_resolution,
            dtype="f4",
        )
        okbins = sorted(okbins)
        oklon = self.compbins["lon"][okbins]
        for val in ("weighted_values", "weighted_squares", "counts"):
            full = np.zeros(
                (self.compbins["lat"].size, lonbins.size), dtype="f4"
            )
            base = getattr(self, val)
            full[:, np.isin(lonbins, oklon)] = base[:, okbins]
            for w in outofbound_recs:
                full[:, np.isin(lonbins, w["lon"])] += base[:, w["outbins"]]
            setattr(self, val, full)
        self.fullbins = {"lon": lonbins, "lat": self.compbins["lat"].copy()}

    @only_when(attr="_running")
    def _update_runmessage(self):
        if (n_ready := len(self.done)) < 1:
            return
        s_per_tile = self.watch.total / n_ready
        eta = (len(self.argrecs) - n_ready) * s_per_tile / 3600
        with self._setmode():
            self.runmessage = (
                f"{n_ready}/{len(self.argrecs)} tiles "
                f"({len(self.failures)} failed) | "
                f"({round(s_per_tile, 2)} s/tile; "
                f"{round(self.watch.total / 3600, 2)} h elapsed; "
                f"est. {round(eta, 2)} h remaining)"
            )

    @only_when(attr="_running")
    def _update_from_tile(self, i, result):
        tbins = result.pop("tbins")
        lats, lons = tbins["lat_bins"], tbins["lon_bins"]
        for a in ("weighted_values", "weighted_squares", "counts"):
            getattr(self, a)[slice(*lats), slice(*lons)] += result.pop(a)
        if self.extended_info is True:
            self.scounts[slice(*lats), slice(*lons)] += result.pop('scounts')
            # for set1, set2 in zip(
            #     self.orbs[slice(*lats), slice(*lons)].ravel(),
            #     result.pop('orbs').ravel()
            # ):
            #     set1.update(set2)
        self.result_meta[i] = result

    @only_when(attr="_running")
    def _pollfetch(self):
        for k in self.ready_results():
            try:
                v = self.pool.results.pop(k).get()
                self._update_from_tile(k, v)
            except KeyboardInterrupt:
                raise
            except Exception as ex:
                self.failures[k] = ex
            self.done.add(k)

    @only_when(attr="_running")
    def _execute_serial(self, verbose):
        with self._setmode():
            self.runmessage = "no tiles complete"
        for i, rec in enumerate(self.argrecs):
            try:
                self._update_from_tile(i, process_tile(**rec['kwargs']))
            except KeyboardInterrupt:
                raise
            except Exception as ex:
                self.failures[i] = ex
            self.done.add(i)
            self.watch.click()
            self._update_runmessage()
            if verbose is True:
                print_inline(self.runmessage)
        with self._setmode():
            self._done = True
            self._finalize_run()

    @only_when(attr="_running")
    def _execute(self, verbose):
        with self._setmode():
            self.runmessage = "no tiles complete"
        try:
            while not len(self.done) == len(self.argrecs):
                self._pollfetch()
                self.watch.click()
                self._update_runmessage()
                if verbose is True:
                    print_inline(self.runmessage)
                if len(self.done) != len(self.argrecs):
                    time.sleep(0.25)
            with self._setmode():
                self._done = True
                self._finalize_run()
        finally:
            self.pool.terminate()

    @only_when(attr="_done")
    def _finalize_run(self):
        self._rectify_coordinates()
        if self.mapped_gain_cutoff is not None:
            cut = self.counts < self.mean_case_count_threshold()
            calc_counts = self.counts.copy()
            calc_counts[cut] = 0
            if self.extended_info is True:
                self.scounts[cut] = 0
        else:
            calc_counts = self.counts
        self.temp = self.weighted_values / calc_counts
        self.std = biased_std(
            self.weighted_values, self.weighted_squares, calc_counts
        )
        self.valid, self.vslice, self.tempcut = cut_valid(self.temp)
        self._make_backplanes()

    @only_when(attr="_prepped")
    def run(self, n_threads=None, verbose=True, block=True):
        if self._armed is False or n_threads != self.pool.threads:
            self.arm(n_threads)
        with self._setmode():
            self._running, self._armed = True, False
            self.done = set()
        self.watch.click()
        with self._setmode():
            if self.pool.threads is not None:
                self.pool.map(process_tile, self.argrecs)
                if block is False:
                    return ThreadPoolExecutor(1).submit(self._execute, False)
                return self._execute(verbose)
            if block is False:
                warnings.warn("n_threads=None; block=False ignored")
            return self._execute_serial(verbose)

    def print_message(self):
        if hasattr(self, "runmessage"):
            print_inline(self.runmessage)

    @only_when(attr="_prepped")
    def arm(self, n_threads: Optional[int] = None):
        with self._setmode():
            self._armed = True
            self.pool = MaybePool(n_threads)
            self.watch = Stopwatch()
            self.weighted_values = np.zeros(
                (self.compbins["lat"].size, self.compbins["lon"].size),
                dtype="f4",
            )
            self.counts = self.weighted_values.copy()
            self.weighted_squares = self.weighted_values.copy()
            if self.extended_info is True:
                self.tilecounts = np.zeros(
                    self.weighted_values.shape, dtype=np.uint16
                )  # EXTENDED
                self.scounts = np.zeros(
                    self.weighted_values.shape, dtype=np.uint32
                )  # EXTENDED
                # self.orbs = np.full(
                #     self.weighted_values.shape, None, dtype="object"
                # )  # E XTENDED
                # for i, j in np.ndindex(self.orbs.shape):  # EXTENDED
                #     self.orbs[i, j] = set()
            self.failures, self.result_meta = {}, {}

    @only_when(attr="_prepped")
    def project_power_pattern(self, d, lat, lon, bins=None):
        bins = bins if bins is not None else self.compbins
        a2l_constants = get(list(A2L_CONSTANTS), self.constants)
        plat, plon = self.a2l(d, lat, lon, *a2l_constants)
        # noinspection PyTypeChecker,PyUnresolvedReferences
        return binned_statistic_2d(
            plat,
            plon,
            self.constants['WEIGHTS'].ravel(),
            bins=(bins['lat'].size, bins['lon'].size),
            range=(
                (bins['lat'].min(), bins['lat'].max()),
                (bins['lon'].min(), bins['lon'].max()),
            )
        ).statistic

    # noinspection PyUnboundLocalVariable
    @only_when(attr="_prepped")
    def mean_case_count_threshold(self):
        """
        this is a heuristic. in fact there's a complex
        relationship between lat and height. but whatever.
        """
        if self.mapped_gain_cutoff is None:
            return None
        if len(self.table) == 0:
            raise ValueError("Empty table, cannot calculate mean cutoff.")
        if not hasattr(self, "tilespecs") or len(self.tilespecs) == 0:
            raise ValueError("No tilespecs available.")
        latsort = sorted(self.tilespecs, key=lambda t: max(map(abs, t['lat'])))
        cases = []
        for direction in(iter, reversed):
            for tile in direction(latsort):
                tileframe = load_tileframe(
                    tile,
                    tablepath=self.fast_tablepath,
                    presliced=True,
                    channel=self.channel
                )
                if len(tileframe) > 0:
                    break
            worst = tileframe.loc[
                tileframe['lat'] == tileframe['lat'].max()
            ].iloc[[0]]
            a2l_constants = get(list(A2L_CONSTANTS), self.constants)
            tbins = get_tile_bins(
                self.compbins, tile, worst['d'], self.a2l, a2l_constants
            )
            oneshot = integrate_tile(
                worst,
                tbins,
                self.a2l,
                self.sigma,
                self.constants['WEIGHTS'],
                a2l_constants
            )
            power = self.project_power_pattern(
                *worst.iloc[0][['d', 'lat', 'lon']], tbins
            )
            cases.append(
                oneshot['counts'][
                    np.nonzero(power >= db2power(self.mapped_gain_cutoff))
                ].min()
            )
        return np.mean(cases)

    def __str__(self):
        selfstr = f"DeconvManager ({self.state})\n"
        if hasattr(self, "argrecs"):
            selfstr += f"  {len(self.argrecs)} total tiles\n"
        if hasattr(self, "done"):
            selfstr += f"  {len(self.done)} tiles complete\n"
        if hasattr(self, "failures"):
            if len(self.failures) > 0:
                selfstr += f"  ({len(self.failures)} failed)\n"
        if hasattr(self, "pool"):
            if self.pool.threads is not None:
                selfstr += f"  {self.pool.threads} threads\n"
            else:
                selfstr += f"  serial\n"
        for k, v in self.settings.items():
            rep = str(v)
            if len(rep) > 30:
                rep = f"{rep[:10]}...{rep[-10:]}"
            selfstr += f"  {k}: {rep}\n"
        return selfstr

    def __repr__(self):
        return self.__str__()

    @property
    def settings(self):
        return {
            attr: getattr(self, attr)
            for attr in getfullargspec(DeconvManager.__init__).kwonlyargs
        }

    @property
    def state(self):
        for state in reversed(self.states):
            if getattr(self, state) is True:
                return state.strip("_")
        return "not prepped"

    def copy(self, **kwargs):
        return DeconvManager(**(self.settings | kwargs))

    states = ("_prepped", "_armed", "_running", "_done")
    runattrs = (
        "pool",
        "watch",
        "weighted_values",
        "counts",
        "weighted_squares",
        "failures",
        "result_meta",
        "runmessage",
        "temp",
        "tempcut",
        "fullbins",
    )
    prepattrs = (
        "chunks",
        "tilespecs",
        "argrecs",
        "a2l",
        "constants",
        "tpath",
        "table",
        "max_abslat",
        "max_orbalt",
        "spill",
        "argrecs",
        "max_orbalt",
        "compbins",
        "nb_cachename",
    )

    _armed = False
    _prepped = False
    _readied = False
    _running = False
    _done = False
    a2l: Callable
    abslat_cutoff: float
    backplanes: dict[str, np.ndarray]
    constants: dict
    chunks: tuple[tuple[float, float], ...]
    done: set
    failures: dict
    compbins: dict[str, np.ndarray]
    fullbins: dict[str, np.ndarray]
    integration_constants = (
        "WEIGHTS",
        "SCALED_RAD",
        "Y_OVER_RAD",
        "X_OVER_RAD",
    )
    max_abslat: np.float32
    max_orbalt: np.float32
    weighted_values: np.ndarray
    weighted_squares: np.ndarray
    counts: np.ndarray
    nb_cachename: str
    std: np.ndarray
    pool: MaybePool
    result_meta: dict
    runmessage: str
    spill: dict[str, np.float32]
    temp: np.ndarray
    tempcut: np.ndarray
    tilespecs: list[dict]
    fast_tablepath: Path
    table: pa.Table
    valid: tuple[np.ndarray, np.ndarray] = None
    vslice: tuple[slice, slice]
    # extended attributes
    tilecounts: np.ndarray
    scounts: np.ndarray
    orbs: np.ndarray
    watch: Stopwatch()
