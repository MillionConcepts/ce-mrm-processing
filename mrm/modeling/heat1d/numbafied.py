import math

import numba as nb
import numpy as np


def heat_capacity_scalar(polynomial: np.ndarray, temperature: float):
    output = 0
    for i, term in enumerate(polynomial):
        output += temperature ** (len(polynomial) - i - 1) * term
    return output


def heat_capacity(polynomial: np.ndarray, temperature: np.ndarray):
    output = np.zeros(temperature.shape)
    for i, term in enumerate(polynomial):
        output += temperature ** (len(polynomial) - i - 1) * term
    return output


def thermcond(temp: np.ndarray, kc: np.ndarray, r350: float):
    return kc * (1 + r350 * temp**3)


def update_thermals(
    temp: np.ndarray, kc: np.ndarray, r350: float, polynomial: np.ndarray
):
    cp = heat_capacity(polynomial, temp)
    k = thermcond(temp, kc, r350)
    # k = np.full_like(temp, 0.0025)
    return cp, k


def surftemp(
    temp: nb.float64[:],
    dtsurf: float,
    emissivity: nb.float64,
    sigma: nb.float64,
    dz: np.ndarray,
    kc: np.ndarray,
    r350: float,
    qs: float,
) -> nb.float64:
    ts = temp[0]
    delta_t = ts
    # NOTE: C scales the cutoff criterion like this
    cutoff = dtsurf * ts
    # NOTE: commented-out code restates the expressions in the same form
    #  they appear in the C. results are identical.
    twodz = 2 * dz[0]
    layer_0_b = kc[0] * r350
    while np.abs(delta_t) > cutoff:
        dtdz = (-3 * ts + 4 * temp[1] - temp[2]) / twodz
        rad = emissivity * sigma * ts**4
        k = kc[0] + layer_0_b * ts**3  # thermCond
        f1 = rad - qs - k * dtdz
        f2 = (
            4 * rad / ts
            + 3 * (kc[0] + layer_0_b * ts**3) / twodz
            - 3 * layer_0_b * ts**2 * dtdz
        )
        delta_t = -f1 / f2
        ts = ts + delta_t
    return ts


def update_t(
    temp: np.ndarray,
    k: np.ndarray,
    cp: np.ndarray,
    surfflux: float,
    kc: np.ndarray,
    g1,
    g2,
    rho: np.ndarray,
    dt: float,
    dz: np.ndarray,
    r350: float,
    qb: float,
    emissivity: float,
    sigma: float,
    dtsurf: float,
) -> np.ndarray:
    c1, c2, c3 = make_propagation_parameters(g1, g2, k, temp)
    for i in range(1, len(temp) - 2):
        temp[i] += (
            dt
            * (c1[i - 1] * temp[i - 1] + c2[i - 1] * temp[i] + c3[i - 1] * temp[i + 1])
            / (rho[i] * cp[i])
        )
    temp[0] = surftemp(temp, dtsurf, emissivity, sigma, dz, kc, r350, surfflux)
    kp = (k[-2] + k[-3]) / 2
    temp[-2] = temp[-3] + (qb / kp) * dz[-1]
    return temp


def make_propagation_parameters(g1, g2, k, temp):
    c1, c2, c3 = (
        np.zeros(len(temp) - 2),
        np.zeros(len(temp) - 2),
        np.zeros(len(temp) - 2),
    )
    for i in range(0, len(temp) - 2):
        c1[i] = g1[i] * k[i]
        c3[i] = g2[i] * k[i + 1]
        c2[i] = -1 * (c1[i] + c3[i])
    return c1, c2, c3


def surf_flux(
    t, dec, r, day, lat, albedo, albedo_coef_0, albedo_coef_1, sabs, rau
):
    h = 2 * np.pi * t / day  # hour angle
    # science!
    cos_i = np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(h)
    # clipping function -- set incidence to 0 when sun below horizon
    cos_i = 0.5 * (cos_i + np.abs(cos_i))
    i = np.arccos(cos_i)
    albedo_model = (
        albedo
        + albedo_coef_0 * (i / np.pi * 4) ** 3
        + albedo_coef_1 * (i / np.pi * 2) ** 8
    )
    f = (1 - albedo_model) / (1 - albedo)
    # NOTE: C code adds another term to cos_i but it appearts to always be 0
    #  unless slope is set nonzero, and by default it is not
    # NOTE: no explicit r / rau term in C
    return f * sabs * (r / rau) ** -2 * cos_i


def orbit_params(
    nu, sin_obliq, useful_parameter, another_useful_parameter, ecc, lp, dt
):
    # Distance to Sun
    r = useful_parameter / (1 + ecc * math.cos(nu))
    # Solar declination
    dec = math.asin(sin_obliq * math.sin(nu + lp))
    # r ** -2 * another_useful_parameter is angular velocity
    return r, dec, nu + r**-2 * dt * another_useful_parameter


# this ugly-looking pair of functions is intended to spare us the very, very
# intense overhead of calling advance() over and over from the Python
# interpreter.
def advance_to_equilibrium(
    albedo,
    albedo_coef_1,
    albedo_coef_2,
    another_useful_parameter,
    cp,
    cpcoeff,
    day,
    dt,
    dtsurf,
    dz,
    eccentricity,
    emissivity,
    equiltime,
    g1,
    g2,
    k,
    kc,
    lat,
    lp,
    nu,
    qb,
    r350,
    rau,
    rho,
    sabs,
    sigma,
    sin_obliq,
    t,
    temp,
    useful_parameter,
):
    while t < equiltime:
        temp, k, cp, nu = advance(
            albedo,
            albedo_coef_1,
            albedo_coef_2,
            another_useful_parameter,
            cp,
            cpcoeff,
            day,
            dt,
            dtsurf,
            dz,
            eccentricity,
            emissivity,
            equiltime,
            g1,
            g2,
            k,
            kc,
            lat,
            lp,
            nu,
            qb,
            r350,
            rau,
            rho,
            sabs,
            sigma,
            sin_obliq,
            t,
            temp,
            useful_parameter,
        )
        t += dt
    return temp, k, cp, nu, t


def advance(
    albedo,
    albedo_coef_1,
    albedo_coef_2,
    another_useful_parameter,
    cp,
    cpcoeff,
    day,
    dt,
    dtsurf,
    dz,
    eccentricity,
    emissivity,
    _equiltime,
    g1,
    g2,
    k,
    kc,
    lat,
    lp,
    nu,
    qb,
    r350,
    rau,
    rho,
    sabs,
    sigma,
    sin_obliq,
    t,
    temp,
    useful_parameter,
):
    r, dec, nu = orbit_params(
        nu,
        sin_obliq,
        useful_parameter,
        another_useful_parameter,
        eccentricity,
        lp,
        dt,
    )
    flux = surf_flux(
        t,
        dec,
        r,
        day,
        lat,
        albedo,
        albedo_coef_1,
        albedo_coef_2,
        sabs,
        rau,
    )
    temp = update_t(
        temp,
        k,
        cp,
        flux,
        kc,
        g1,
        g2,
        rho,
        dt,
        dz,
        r350,
        qb,
        emissivity,
        sigma,
        dtsurf,
    )
    cp, k = update_thermals(temp, kc, r350, cpcoeff)
    return temp, k, cp, nu


nb.jit_module(cache=True, nopython=True)
