from types import MappingProxyType as MPt

import numba as nb
import numpy as np
from quickseries import quickseries

from mrm.shared.constants import HALFPI, LUNAR_RADIUS_KM, PI, RADIANS, TWOPI


# NOTE: if we need to get more specific, we can compute the ray-sphere
# intersection for each grid cell, but this will be much slower at 
# _several_ steps and probably add very little to the result. 

# also, this all gets weird at the poles -- we'll have to use a
# different projection if we want polar maps.

def km_offset_footprint(
    rix_y: nb.float32[:],  # centered pixel index y-axis raveled grid
    rix_x: nb.float32[:],  # centered pixel index x-axis raveled grid
    orbalt: nb.float32,  # orbital altitude (km)
    npix: nb.float32,  # number of pixels on each grid axis
    altaz_extent_rad: nb.float32  # angular half-extent of viewport in radians
) -> tuple[nb.float32[:], nb.float32[:]]:
    """
    convert beam alt/az indices to offset from boresight on surface. 
    this ignores the effect of curvature due to the small extent of the beam;
    i.e., it treats the beam as a simple orthographic projection onto the 
    lunar surface.
    
    if we need to get more specific, we can compute the ray-sphere 
    intersection for each grid cell, but this will be much slower at 
    _several_ steps and probably add very little to the result. note that 
    this _will_ get a little bit weird at the poles, but so 
    will most things.
    """
    # just treat this like a triangle:
    # tangent of beam_angle = half of total surface extent divided by orbital
    # altitude;
    # then surface extent = tan(view_angle / 2) * orbalt * 2
    km_per_pixel = np.tan(altaz_extent_rad) * orbalt * 2 / npix
    # offset at center of pixel from boresight, in km
    return rix_y * km_per_pixel, rix_x * km_per_pixel


# components for piecewise-approximated arcsine function. note that this is
# actually slower _on its own_ in many cases than np.arcsin. however,
# it is faster inside other numbafied function, presumably because it can
# reuse cached pointers.
ASIN_KWARGS = MPt(
    {'jit': True, 'fitres': 40, 'func': 'asin(x)', 'precision': 32}
)
asin1 = quickseries(nterms=8, bounds=(0, 0.5), **ASIN_KWARGS)
asin2 = quickseries(nterms=7, bounds=(0.5, 0.75), **ASIN_KWARGS)


@nb.njit
def qasin(x: nb.float32[:]) -> nb.float32[:]:
    output = np.empty_like(x)
    for i in range(len(x)):
        s = np.sign(x[i])
        el = s * x[i]
        if el < np.float32(0.5):
            output[i] = s * asin1(el)
        elif el < np.float32(0.75):
            output[i] = s * asin2(el)
        else:
            output[i] = np.arcsin(s * el)
    return output


# numbafied quickseries versions of simple trig functions.
# simple trig functions are so optimized in standard libraries that most of
# these offer only small improvements over stock numpy versions if not jitted
# (1.3-3x). However, unlike stock numpy simple trig functions, they _can_ be
# usefully jitted! Don't use these as components of more complicated
# univariate functions involving multiple trig functions
# (e.g. cos(x) + sin(x) * cos(x**2)). You will get much better performance if
# you run quickseries on the more complicated function as a whole.

qatan = quickseries(
    "atan(x)", nterms=14, bounds=(-1.5, 1.5), jit=True, precision=32
)
"""~19x improvement over 32-bit np.arctan."""


qtan = quickseries(
    "tan(x)", nterms=10, bounds=(-1.5, 1.5), jit=True, precision=32
)
"""~40x improvement over 32-bit np.tan."""

qsin = quickseries(
    "sin(x)", bounds=(0, 1), jit=True, nterms=7, precision=32
)
"""
~5x improvement over 32-bit np.sin. Works well from -1 to 1, and not horribly 
to about +/- 1.4.
"""

# same comments as qsin, except bounds of okness are more like +/- 1.2
qcos = quickseries(
    "cos(x)", bounds=(0, 1), jit=True, nterms=6, precision=32
)
"""
~5x improvement over 32-bit np.cos. Works well from -1 to 1, and not horribly 
to about +/- 1.2.
"""


@nb.njit
def qatan2(y0: nb.float32[:], x0: nb.float32[:]) -> nb.float32[:]:
    """numbafied atan2 reimplementation. ~9x speedup over np.arctan2."""
    base = qatan(y0/x0)
    # assign to quadrants
    for i in range(len(y0)):
        xsign = np.sign(x0[i])
        if xsign == 1:
            continue
        ysign = np.sign(y0[i])
        if xsign == -1:
            if ysign == -1:
                base[i] -= PI
            else:
                base[i] += PI
        # and when x is 0, pick the special value to assign
        elif ysign == 0:
            base[i] = 0
        else:
            base[i] = HALFPI * ysign
    return base


@nb.njit
def azel_to_latlon(
    rix_y: nb.float32[:],
    rix_x: nb.float32[:],
    rix_rad: nb.float32[:],
    orbalt: nb.float32,
    npix: nb.float32,
    altaz_extent_rad: nb.float32,
    lat_0: nb.float32,
    lon_0: nb.float32
) -> tuple[nb.float32[:], nb.float32[:]]:
    """
    Accelerated version of reverse orthographic projection function. Projects
    instrument footprint, measured in az/el relative to boresight, to lat/lon
    on planetary surface. Provides ~2x speedup over unaccelerated version,
    slightly more at lower latitudes and slightly less at higher.

    Approximates the surface as locally flat and assumes the instrument is
    looking perfectly nadir. This allows us to treat the offset from center
    as a right triangle, so linear surface extent in rectangular coordinates
    is simply tan(angular offset from boresight / 2) * orbital altitude * 2.

    Then we assume that the target body is perfectly spherical and transform
    these rectangular coordinates into lat/lon using standard spherical trig
    formulae.
    """

    # kilometers per pixel
    kpp = np.tan(altaz_extent_rad) * orbalt * np.float32(2) / npix
    c = qasin(rix_rad * kpp / LUNAR_RADIUS_KM)
    cos_c, sin_c = qcos(c), qsin(c)
    lat = qasin(
        cos_c * np.sin(lat_0)
        + (kpp*rix_y * sin_c * np.cos(lat_0))
        / (rix_rad * kpp + np.float32(1e-8))
    )
    lon = lon_0 + qatan2(
        kpp * rix_x * sin_c,
        rix_rad * kpp * cos_c * np.cos(lat_0)
        - kpp * rix_y * sin_c * np.sin(lat_0)
    )
    for i in range(len(lon)):
        if lon[i] >= TWOPI:
            lon[i] = (lon[i] / RADIANS) - np.float32(360)
        elif lon[i] < 0:
            lon[i] = (lon[i] / RADIANS) + np.float32(360)
        else:
            lon[i] = lon[i] / RADIANS
    return lat / RADIANS, lon


def maxdiff(ref, test, *args):
    """utility arithmetic error-checking function"""
    refval = ref(*args)
    diff = np.abs(refval - test(*args))
    rdiff_raw = diff / refval
    # it's not pathological for refval to be 0, etc.
    rdiff = np.abs(rdiff_raw[np.isfinite(rdiff_raw)])
    return diff.max(), rdiff.max() * 100


def coef_det(fit_curve: np.ndarray, dependent: np.ndarray) -> float:
    """coefficient of determination"""
    # sample variance
    ss_t = sum((dependent - dependent.mean()) ** 2)
    # sum of squares of residuals
    ss_r = sum((dependent - fit_curve) ** 2)
    det = 1 - ss_r / ss_t
    return det
