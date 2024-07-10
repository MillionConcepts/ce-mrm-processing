"""Physical, mathematical, and engineering constants."""

from types import MappingProxyType as MPt

from astropy import constants
import numpy as np

GM = 3.96423e-14  # G*Msun [AU**3/s**2]


# 32-bit versions of some constants to avoid upcasting. The silly ones like PI
# are for readability in numbafied functions.

C = np.float32(constants.c.value)
PLANCK = np.float32(constants.h.value)
BOLTZMANN = np.float32(constants.k_B.value)
SIGMA_SB = np.float32(constants.sigma_sb.value)
LUNAR_RADIUS = np.float32(1737400)
LUNAR_RADIUS_KM = np.float32(LUNAR_RADIUS / 1000)
RADIANS = np.float32(np.pi / 180)
PI = np.float32(np.pi)
TWOPI = np.float32(2 * PI)
HALFPI = np.float32(np.pi / 2)

ANTENNA_PARAMS = {
    # 3 GHz
    1: {
        "mainbeamwidth": 32.2,
        "sigma": 7.6,
        "m_power": -16.5,
        "sidepeak": -14.8,
        "s_by_m": 0.53,
    },
    # 7.8 GHz
    2: {
        "mainbeamwidth": 24,
        "sigma": 4.5,
        "m_power": -36,
        "sidepeak": -16.7,
        "s_by_m": 0.67,
    },
    # 19.35 GHz
    3: {
        "mainbeamwidth": 24,
        "sigma": 5.5,
        "m_power": -18,
        "sidepeak": -17,
        "s_by_m": 0.2,
    },
    # 37 GHz
    4: {
        "mainbeamwidth": 26,
        "sigma": 5.3,
        "m_power": -30,
        "sidepeak": -20,
        "s_by_m": 0.23,
    },
}
"""Canonical antenna pattern model parameters by MRM channel."""

CHANNEL_FREQ = MPt({1: 3, 2: 7.8, 3: 19.35, 4: 37})
"""Nominal center frequency, in GHz, for each MRM channel."""
