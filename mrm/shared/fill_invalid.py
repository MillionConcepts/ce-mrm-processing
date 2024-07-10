import numpy as np
from scipy.interpolate import griddata


def reg_interp(df, ax='et', var='t1', step=1):
    x = np.arange(df[ax].min(), df[ax].max() + 1, step)
    return x, np.interp(x, df[ax], df[var])


def fill_invalid(arr, method='linear', inplace=False):
    val = np.isfinite(arr)
    if val.all():
        return arr
    x, y = np.indices(arr.shape)
    vx, vy, nvx, nvy = x[val], y[val], x[~val], y[~val]
    del val
    canvas = arr if inplace is True else arr.copy()
    canvas[nvx, nvy] = griddata(
        (vx, vy),
        canvas[vx, vy],
        (nvx, nvy),
        method=method
    )
    return canvas
