from itertools import product

import numba as nb
import numpy as np
import pandas as pd
import pyarrow as pa
from pyarrow import compute as pc

from mrm.shared.constants import C, PLANCK, BOLTZMANN


def std_encode(arr, keytype=pa.uint16(), valtype=None, length=None):
    """Inline dictionary encoding for recursively-assembled pyarrow tables."""
    if not isinstance(arr, (pa.Array, list, tuple, np.ndarray, pd.Series)):
        arr = pa.array([arr for _ in range(length)])
    elif not isinstance(arr, pa.Array):
        arr = pa.array(arr)
    if valtype is None:
        valtype = arr.type
    return pc.dictionary_encode(arr).cast(pa.dictionary(keytype, valtype))


@nb.jit(cache=True, nopython=True)
def simulate_brightness_temperature(
    temp: nb.float32[:],
    dz: nb.float32[:],
    rho: nb.float32[:],
    frequency_ghz: nb.float32(),
    tio2_percentage: nb.float32(),
    a: nb.float32(),
    b: nb.float32(),
    c: nb.float32(),
) -> tuple[nb.float32[:], nb.float32[:]]:
    """
    calculate brightness temperature in a non-scattering medium
    based on physical temperature, layer thickness, layer density,
    observation frequency, and material parameters.

    temp.shape[0] must be equal to dz.size & rho.size. Additional
    columns of temp are assumed to share dz / rho values, and may
    represent temperature readings at different times, etc. --
    per-column calculations are entirely independent of one another
    """
    freq = frequency_ghz * 10 ** 9
    loss_tangent = 10 ** (a * tio2_percentage + b * rho / 1000 - c)
    permittivity = 1.92 ** (rho / 1000)
    attenuation_coefficient = (
            2 * np.pi * freq / C * np.sqrt(permittivity) * loss_tangent
    )
    attenuation = np.exp(-rho * attenuation_coefficient)
    # intra-layer reflection
    lower = np.sqrt(permittivity[1:])
    upper = np.sqrt(permittivity[:-1])
    reflection = ((lower - upper) / (lower + upper)) ** 2
    reflection = np.append(reflection, 0)
    surface_reflectivity = ((upper[0] - 1) / (upper[0] + 1)) ** 2
    # multiple reflection factor
    surface_multiple = (
        1 - surface_reflectivity * reflection[0] * attenuation[0] ** 2
    )
    lower_multiple = (
        1 - reflection[0:-1] * reflection[1:] * attenuation[1:] ** 2
    )
    multiple = np.append(np.array([surface_multiple]), lower_multiple)
    # weights for contribution of each layer to brightness temperature
    weights = np.zeros(temp.shape[0])
    for i in range(temp.shape[0]):
        weights[i] = (
            (1 - attenuation[i])
            * (1 + reflection[i] * attenuation[i])
            * (1 - surface_reflectivity)
            * (1 - reflection[:i]).prod()
            * attenuation[:i].prod()
            / multiple[:i].prod()
        )
    # per-layer brightness temperatures
    ib_temp = np.zeros_like(temp)
    for i in range(temp.shape[1]):
        ib_temp[:, i] = (
                1 / (np.exp(
            PLANCK * freq / BOLTZMANN / temp[:, i]) - 1) * weights
        )
    # base brightness temperature
    tb = PLANCK * freq / BOLTZMANN / np.log(1 / np.sum(ib_temp, axis=0) + 1)
    # TB contribution from layers 'below the bottom' (??)
    tr = (
        temp[-1, :]
        + (temp[-1, :] - temp[-2, :])
        * 2
        / (dz[-2] + dz[-1])
        / attenuation_coefficient[-1]
    )
    # The TB from the layers 'below the bottom' after attenuation
    t_bottom = (
        tr
        * (1 - surface_reflectivity)
        * (1 - reflection).prod()
        * attenuation.prod()
        / multiple.prod()
    )
    return tb + t_bottom, weights


def simulate_brightnesses(
    lat,
    values,
    ltst,
    frequencies,
    tio2_percentage,
    a,
    b,
    c_values
):
    print(f"running {round(float(lat), 1)}")
    tb_results = []
    for alb_h, values in values.groupby(['albedo', 'h']):
        alb, h = alb_h
        temps = []
        for lt, timegroup in values.groupby('ltst'):
            temps.append(timegroup['temp'])
        temp = np.vstack(temps).T
        # noinspection PyUnboundLocalVariable
        dz, rho = timegroup['dz'].values, timegroup['rho'].values
        del temps
        for freq, c in product(frequencies, c_values):
            # these values should be presorted by the assemble_heat1d_outputs
            # pipeline -- could add a check if there is concern
            tb, _ = simulate_brightness_temperature(
                temp, dz, rho, freq, tio2_percentage, a, b, c
            )
            arrays = {
                'freq': std_encode(freq, pa.uint8(), pa.float32(), len(ltst)),
                'c': std_encode(c, pa.uint8(), length=len(ltst)),
                'lat': std_encode(lat, pa.uint8(), length=len(ltst)),
                'albedo': std_encode(alb, pa.uint8(), length=len(ltst)),
                'h': std_encode(h, pa.uint8(), length=len(ltst)),
                'ltst': std_encode(ltst),
                'tb': pa.array(tb).cast(pa.float32())
            }
            # noinspection PyTypeChecker
            table = pa.Table.from_arrays(
                list(arrays.values()), list(arrays.keys())
            )
            tb_results.append(table)
    concat = pa.concat_tables(tb_results).unify_dictionaries().combine_chunks()
    print(f"completed {round(float(lat), 1)}")
    return concat
