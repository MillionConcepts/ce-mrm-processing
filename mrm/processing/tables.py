from collections import defaultdict
from itertools import repeat
from typing import Sequence, Union

import astropy.time as at
from dustgoggles.structures import MaybePool
from hostess.monitors import Stopwatch
from hostess.profilers import Profiler
import numpy as np
from numpy.linalg import norm
import pandas as pd
import pyarrow as pa
from pyarrow import compute as pac
from quickseries.simplefit import fit
import spiceypy as spice

from moonbow.expressions import rowmul
from scipy.stats import median_abs_deviation


PROFILER = Profiler({'time': Stopwatch()})
ESETTINGS = {}


def init_spice_kernel_pool(metakernel):
    from lhorizon.kernels import load_metakernel
    load_metakernel(metakernel)


def utc2et_chunk(utc_times, metakernel):
    init_spice_kernel_pool(metakernel)
    return np.array(list(map(spice.utc2et, utc_times)))


def utc2et_map(chunks, threads):
    pool = MaybePool(threads)
    try:
        pool.map(
            utc2et_chunk, [
                {'args': (c, ESETTINGS['metakernel'])} for c in chunks
            ]
        )
        pool.close()
        pool.join()
        results = pool.get()
        return np.concatenate(list(results.values()))
    finally:
        pool.terminate()


def latsrf_chunk(et_times, longitudes, latitudes, metakernel):
    points = []
    longitudes, latitudes = map(np.radians, (longitudes, latitudes))
    init_spice_kernel_pool(metakernel)
    for et, lon, lat in zip(et_times, longitudes, latitudes):
        srf = spice.latsrf(
            "ELLIPSOID", "MOON", et, "MOON_ME", [[lon, lat]]
        )
        points.append(srf)
    return np.vstack(points)


def latsrf_map(et, lon, lat, threads):
    pool = MaybePool(threads)
    nchunks = 1 if threads is None else threads
    try:
        et_chunks = np.array_split(et, nchunks)
        lon_chunks = np.array_split(lon.values, nchunks)
        lat_chunks = np.array_split(lat.values, nchunks)
        pool.map(
            latsrf_chunk, 
            [
                {'args': (et, lon, lat, ESETTINGS['metakernel'])}
                for et, lon, lat in zip(et_chunks, lon_chunks, lat_chunks)
            ]
        )
        del et_chunks, lon_chunks, lat_chunks
        pool.close()
        pool.join()
        return np.vstack(list(pool.get().values()))
    finally:
        pool.terminate()


def jd_tdb_chunk(utc_times):
    at_time = at.Time(list(utc_times), scale='utc')
    jd = np.array(at_time.jd)
    tdb = pd.Series(at_time.tdb.datetime)
    return jd, tdb


def jd_tdb_map(chunks, threads):
    pool = MaybePool(threads)
    try:
        pool.map(jd_tdb_chunk, [{'args': (c,)} for c in chunks])
        pool.close()
        pool.join()
        results = pool.get()
        jd = np.concatenate([r[0] for r in results.values()])
        tdb = pd.concat(
            [r[1] for r in results.values()]
        ).reset_index(drop=True)
        return jd, tdb
    finally:
        pool.terminate()


def spkcpo_chunk(    
    et_times: Sequence[float],
    observer_positions: Sequence[Sequence[float]] | Sequence[float],
    target_name: Union[int, str],
    center_name: Union[int, str],
    outframe: str,
    inframe: str,
    metakernel: str,
    abcorr: str = "LT+S",
    method: str = "OBSERVER",
):
    """
    compute positions of target (in an arbitrary frame) relative to observers
    with constant positions (defined in a body-fixed frame) at times et_times.
    """
    target_positions = []
    target, center = map(str, (target_name, center_name))
    if isinstance(observer_positions[0], (int, float)):
        observer_positions = repeat(observer_positions)
    init_spice_kernel_pool(metakernel)
    for moment, pos in zip(et_times, observer_positions):
        state, _lt = spice.spkcpo(
            target, moment, outframe, method, abcorr, pos, center, inframe
        )
        target_positions.append(state[:3])
    return np.vstack(target_positions)


def spkcpo_map(et, positions, threads):
    nchunks = 1 if threads is None else threads
    # permit fixed position
    if isinstance(positions[0], (int, float)):
        pos_chunks = [positions for _ in range(nchunks)]
    else:
        pos_chunks = np.array_split(positions, nchunks)
    et_chunks = np.array_split(et, nchunks)
    settings = {
        'target_name': 10, 
        'center_name': 301, 
        'inframe': 'MOON_ME', 
        'outframe': 'ECLIPJ2000',
        'metakernel': ESETTINGS['metakernel']
    }
    spkcpo_recs = [
        {'kwargs': {'observer_positions': p, 'et_times': e} | settings}
        for e, p in zip(et_chunks, pos_chunks)
    ]
    pool = MaybePool(threads)
    try:
        pool.map(spkcpo_chunk, spkcpo_recs)
        del spkcpo_recs, et_chunks, pos_chunks
        pool.close()
        pool.join()
        # results = pool.get()
        return np.vstack(list(pool.get().values()))
    finally:
        pool.terminate()


def et2lst_chunk(et_times, longitudes, metakernel):
    init_spice_kernel_pool(metakernel)
    times = []
    for et, lon in zip(et_times, np.radians(longitudes)):
        times.append(spice.et2lst(et, 301, lon, 'PLANETOCENTRIC')[3])
    # noinspection PyUnresolvedReferences
    times = pd.to_timedelta(times).total_seconds() / (24 * 60 * 60)
    return times.astype('float32').values


def et2lst_map(et, lon, threads):
    nchunks = 1 if threads is None else threads
    et_chunks = np.array_split(et, nchunks)
    lon_chunks = np.array_split(lon.values, nchunks)
    lst_recs = [
        {'args': chunk, 'kwargs': {'metakernel': ESETTINGS['metakernel']}}
        for chunk in zip(et_chunks, lon_chunks)
    ]
    pool = MaybePool(threads)
    try:
        pool.map(et2lst_chunk, lst_recs)
        del lst_recs, et_chunks, lon_chunks
        pool.close()
        pool.join()
        results = pool.get()
        return np.concatenate(list(results.values()))
    finally:
        pool.terminate()    


def temparray(tab):
    return np.asarray(tab.select([f"t{i + 1}" for i in range(4)]))


def check_obstruction(mxyz, xyz0, r):
    """
    mxyz: vector from source to obstructing body center
    xyz0: vector from source to target
    r: radius of obstructing body

    Call the angle between mxyz and xyz0 theta. The length of the line between
    the contact point and the sub-source point is: tan(theta) * |mxyz|.

    Then, we can check for conditions under which the ray-sphere intersection
    solution becomes undefined in a useful way.

    The limit at theta >= arctan(r / |mxyz|) indicates that the ray is
    off-axis from the obstructing body, so the target is targetable.

    The constraint that |xyz0| - |mxyz| == +/- r is only useful on one side:
    a target can be distant from an obstructing body but still obstructed by
    it. However, if the target is _inside_ the body, we will consider it
    obstructed.

    This function assumes:
    1. The observer is not inside the obstructing body (in which case -r
        would be the appropriate constraint, and it would be impossible to be
        off-axis).
    2. The obstructing body is spherical.
    3. The observer and target are small or distant enough that they can be
        treated as having no extent wrt one another.
    4. The target is small enough relative to the obstructing body that it
        can be treated as having no extent wrt to it.
    """
    # vector magnitudes
    mxyz_mag, xyz0_mag = norm(mxyz, axis=1), norm(mxyz, axis=1)
    magdiff = xyz0_mag - mxyz_mag
    # A spherical body can't block a target if the target is closer to an
    # observer than the body's center, unless the target is inside the body.
    closer, outside = magdiff < 0, magdiff <= r
    # unit vectors
    mxyz_dir, xyz0_dir = rowmul(mxyz, 1 / mxyz_mag), rowmul(xyz0, 1 / xyz0_mag)
    offaxis = np.arccos(
        np.sum(mxyz_dir * xyz0_dir, axis=1)  # i.e, dot product
    ) >= np.arctan(r / mxyz_mag)
    return (closer & outside) | offaxis, closer, outside, offaxis


def lat_tb_model(lat, a, b):
    return a * np.cos(lat) ** b


def fit_latitude(table, channel='t1', cutoff=85, guess=(231, 0.35)):
    lat_fit_table = pd.concat(
        [np.radians(table['lat']), table[channel]], axis=1
    )
    lat_fit_table = lat_fit_table.loc[
        table['lat'].abs() < cutoff
        ].reset_index().copy()
    return fit(
        lat_tb_model,
        [lat_fit_table["lat"]],
        lat_fit_table[channel],
        guess=guess
    )


# A module must insert values into ESETTINGS before calling the following
# functions.


def low_temp(tab):
    return np.any(temparray(tab) < ESETTINGS['min_temp'], axis=1), 0b10


def channel_offset(tab):
    return (
        np.ptp(temparray(tab), axis=1) > ESETTINGS['max_channel_offset'],
        0b100
    )


def bad_orbits(tab):
    return pac.is_in(tab['orbit'], pa.array(ESETTINGS['bad_orbits'])), 0b1000


def small_orbits(tab):
    unflagged = tab.filter(ESETTINGS['flagarr'] == 0)['orbit'].to_numpy()
    orbits, orbfreqs = np.unique(unflagged, return_counts=True)
    return pac.is_in(
        tab['orbit'], pa.array(orbits[orbfreqs < ESETTINGS['min_orbit_size']])
    ), 0b10000


def timedupes(tab):
    utimes, timefreqs = np.unique(tab['time'], return_counts=True)
    return np.isin(tab['time'], utimes[timefreqs > 1]), 0b100000


def exclude_orbits(tab):
    return pac.is_in(
        tab['orbit'], pa.array(ESETTINGS['exclude_orbits'])
    ).to_numpy(), 0b1000000


def latmodel_outliers(tab, channel='t1'):
    orbfits = {}
    unflagged = tab.filter(
        ESETTINGS['flagarr'] == 0
    ).select([channel, 'lat', 'orbit'])
    for orbit in pac.unique(unflagged['orbit']).to_numpy():
        orbslice = unflagged.filter(
            pac.equal(unflagged['orbit'], orbit)
        ).to_pandas()
        if len(orbslice) <= 1:
            continue
        orbslice = orbslice.sort_values(by=channel).reset_index().copy()
        orbslice['lat'] = orbslice['lat'].abs()
        (a, b), _det = fit_latitude(orbslice)
        orbfits[orbit] = {
            'a': a,  # constant term in model
            'b': b,  # cosine exponent term
        }
    orbax = np.array(list(orbfits.keys()))
    aax = np.array([rec['a'] for rec in orbfits.values()])
    bax = np.array([rec['b'] for rec in orbfits.values()])
    a_mad = median_abs_deviation(aax)
    a_med = np.ma.median(aax)
    a_off = np.abs(aax - a_med) > a_mad * ESETTINGS['latmodel_cutoff']
    b_mad = median_abs_deviation(bax)
    b_med = np.ma.median(bax)
    b_off = np.abs(bax - b_med) > b_mad * ESETTINGS['latmodel_cutoff']
    return pac.is_in(
        tab['orbit'], pa.array(orbax[a_off | b_off])
    ).to_numpy(), 0b10000000


def latcut_outliers(tab):
    lat_tiles = np.arange(-90, 90, 1)
    table = tab.filter(
        ESETTINGS['flagarr'] == 0
    ).select(['orbit', 'lat', 't1']).to_pandas()
    latcuts = {}
    for tile in lat_tiles:
        latcuts[tile] = table[
            (table['lat'] > tile) & (table['lat'] < tile + 1)
        ][['t1', 'orbit']]
    cutstats = {
        k: {'mean': v['t1'].mean(), 'std': v['t1'].std()}
        for k, v in latcuts.items()
    }
    suspicious = defaultdict(int)
    for tile, cut in latcuts.items():
        cutpivot = cut.pivot_table(
            values='t1', index='orbit', aggfunc=['mean', 'count']
        )
        cutpivot.columns = cutpivot.columns.droplevel(1)
        for orbit, orbstats in cutpivot.iterrows():
            if orbstats['count'] < 10:
                continue
            if np.abs(
                cutstats[tile]['mean'] - orbstats['mean']
            ) > cutstats[tile]['std'] * 2:
                suspicious[orbit] += 1
    svals = np.array(list(suspicious.values()))
    orbax = np.array(list(suspicious.keys()))
    suspicious_orbits = orbax[
        svals > np.percentile(svals, ESETTINGS['tile_vote_centile'])
    ]
    return pac.is_in(
        tab['orbit'], pa.array(suspicious_orbits)
    ).to_numpy(), 0b100000000

