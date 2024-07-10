"""present for convenience. may not stay."""

from typing import Literal, Optional

from lhorizon.lhorizon_utils import sph2cart, listify
from marslab.imgops.pltutils import set_colorbar_font
import matplotlib as mpl
import matplotlib.font_manager as mplf
import matplotlib.pyplot as plt
from moonbow.expressions import hats
from more_itertools import divide
import numpy as np
import pandas as pd

from mrm.shared.constants import LUNAR_RADIUS_KM


def downslice(table, orbit: int, res: int = 30):
    slic = table.loc[table["orbit"] == orbit]
    return slic[0:-1: int(len(slic) / res)]


def orb3d(
    slic,
    ax=None,
    ref: Literal['bodyfixed', 'inertial'] = 'bodyfixed',
    how: Literal['scatter', 'line', 'surface'] = 'scatter',
    colorfield=None,
    cmap="orange_teal",
    color=(1, 1, 1, 1),
    norm=None,
    scattersize: int = 8
):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
    if ref == 'bodyfixed':
        xyz = sph2cart(slic["lat"], slic["lon"], slic["d"] + LUNAR_RADIUS_KM)
    elif ref == 'inertial':
        xyz = [slic[f'm{a}'] - slic[f'c{a}'] for a in ('x', 'y', 'z')]
    else:
        raise ValueError(f'unknown frame type {ref}')
    if colorfield is not None:
        norm = norm if norm is not None else mpl.colors.Normalize()
        if how in ('line', 'surface'):
            raise ValueError('cannot do a line or surf plot with a colorfield')
        kwargs = {'c': slic[colorfield], 'norm': norm, 'cmap': cmap}
    else:
        kwargs = {'color': color}
    if how == 'scatter':
        render = ax.scatter(*xyz, **kwargs, s=scattersize)
    elif how == 'line':
        render = ax.plot(*xyz, **kwargs)
    elif how == 'surface':
        coords = [
            np.vstack([c, np.flip(c)]).T for c in (a.values for a in xyz)
        ]
        render = ax.plot_surface(*coords, **kwargs)
    else:
        raise ValueError(f'unknown plot type {how}')
    return ax, render


def ce_overview_3d(
    table: pd.DataFrame,
    n_orbits: int = 11,
    moonres: int = 50,
    sunres: int = 10,
    orbres: int = 50,
    moonsize: float = LUNAR_RADIUS_KM * 0.8,
    sunsize: float = LUNAR_RADIUS_KM / 10,
    sundist: float = LUNAR_RADIUS_KM * 2,
    cmap: str = "moonbow",
    orbplot_type: Literal["flower", "static", "none"] = "static",
    plot_moon: bool = True,
    plot_sun: bool = True,
):
    """
    3D plot showing orbital pattern, sun position, etc.
    Args:
        table: slice of mrm_concatenated table as DataFrame
        n_orbits: number of orbits/positions to plot
        cmap: colormap for sun spheres and orbits
        orbres: resolution of orbital path
        moonres: resolution of moon sphere
        sunres: resolution of sun sphere
        moonsize: size of moon sphere
        sunsize: size of sun sphere
        sundist: standoff distance of sun sphere
        orbplot_type: how to plot the orbits: "fixed", "flower", or "none".
            probably turn the moon off for flower.
        plot_moon: plot the moon?
        plot_sun: plot the suns?
    """
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    orbs = [next(d) for d in divide(n_orbits, table["orbit"].unique())]
    for i, orb in enumerate(orbs):
        slic = downslice(table, orb, orbres)
        orbcolor = mpl.colormaps[cmap](i * int(255 / len(orbs)))
        pltkwargs = {'ax': ax, 'color': orbcolor, 'slic': slic}
        if orbplot_type == "static":
            # 'plot' looks better but will _not_ render above the Moon
            orb3d(**pltkwargs, ref="inertial", how="scatter")
        elif orbplot_type == "flower":
            orb3d(**pltkwargs, ref="bodyfixed", how="surface")
        elif orbplot_type is not None:
            raise ValueError(f"don't know plot type {orbplot_type}")
        if plot_sun is True:
            norms = hats(slic[["mx", "my", "mz"]].values)
            sunlat = np.linspace(-90, 90, sunres)
            sunlon = np.linspace(0, 360, sunres)
            sunlat, sunlon = map(np.ravel, np.meshgrid(sunlat, sunlon))
            xsun, ysun, zsun = map(
                lambda coord: coord.reshape(sunres, sunres),
                sph2cart(sunlat, sunlon, sunsize),
            )
            # add standoff
            for compix, comp in zip(range(3), (xsun, ysun, zsun)):
                comp += norms[:, compix].mean() * sundist
            ax.plot_surface(xsun, ysun, zsun, color=orbcolor, shade=True)
    if plot_moon is True:
        lat = np.linspace(-90, 90, moonres)
        lon = np.linspace(0, 360, moonres)
        lat, lon = map(np.ravel, np.meshgrid(lat, lon))
        xsph, ysph, zsph = map(
            lambda coord: coord.reshape(moonres, moonres),
            sph2cart(lat, lon, moonsize),
        )
        moon = ax.plot_surface(
            xsph, ysph, zsph, color=(1, 1, 1, 0.3), shade=False
        )
        moon.set_edgecolor([0, 0, 0, 0.6])
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_axis_off()
    ax.set_facecolor([0, 0, 0, 1])
    fig.tight_layout()
    return fig, ax


def ce_orbit_plot(
    table,
    orbits,
    orbres=50,
    moonsize=LUNAR_RADIUS_KM * 0.8,
    plot_moon=True,
    text_standoff=LUNAR_RADIUS_KM * 1.5,
    lon_interval: int = 30,
    lat_interval: int = 15,
    font="Fira Mono",
    fontsize=12,
    loncolor=(1, 0.8, 0.85, 0.9),
    latcolor=(0.8, 0.85, 1, 0.9),
    colorfield: Optional[str] = "phase",
    cmap: str = "PuOr",
    color_crop: tuple[float, float] = (1, 99),
    scattersize: int = 8
):
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    if colorfield is not None:
        norm = mpl.colors.Normalize(
            *np.percentile(table[colorfield], color_crop)
        )
    else:
        norm, cmap = None, None
    for orbit in listify(orbits):
        slic = downslice(table, orbit, orbres)
        _, mappable = orb3d(
            slic,
            ax,
            colorfield=colorfield,
            cmap=cmap,
            norm=norm,
            scattersize=scattersize
        )
    compass_lat = np.arange(-90, 90 + lat_interval, lat_interval)
    compass_lon = np.arange(-180, 180 + lon_interval, lon_interval)
    latx, laty, latz = sph2cart(compass_lat, 0, text_standoff)
    lonx, lony, lonz = sph2cart(0, compass_lon, text_standoff)
    fconfig = mplf.FontProperties(family=[font], size=fontsize)
    # TODO: would be really nice to have ticks
    if not hasattr(lonz, "len"):
        lonz = [lonz for _ in range(len(lonx))]
    for x0, y0, z0, val in zip(latx, laty, latz, compass_lat):
        ax.text(x0, y0, z0, round(val), c=latcolor, fontproperties=fconfig)
    # don't actually want to draw the 180
    for x0, y0, z0, val in zip(lonx, lony, lonz, compass_lon[:-1]):
        ax.text(x0, y0, z0, round(val), c=loncolor, fontproperties=fconfig)
    # compass lines
    ax.plot(
        *sph2cart(
            np.full(361, 0),
            np.arange(-180, 181, 1),
            np.full(361, text_standoff * 0.85)
        ),
        c=loncolor
    )
    ax.plot(
        *sph2cart(
            np.arange(-90, 91, 1),
            np.full(181, 0),
            np.full(181, text_standoff * 0.9)
        ),
        c=latcolor
    )
    if plot_moon is True:
        # moon sphere object
        mlat, mlon = map(np.ravel, np.meshgrid(compass_lat, compass_lon))
        xsph, ysph, zsph = map(
            lambda coord: coord.reshape(len(compass_lon), len(compass_lat)),
            sph2cart(mlat, mlon, moonsize),
        )
        moon = ax.plot_surface(
            xsph, ysph, zsph, color=(1, 1, 1, 0.1), shade=False
        )
        moon.set_edgecolor([0, 0, 0, 0.6])
    if colorfield is not None:
        # noinspection PyUnboundLocalVariable,PyTypeChecker
        colorbar = fig.colorbar(
            mappable,
            aspect=40,
            shrink=0.8
        )
        set_colorbar_font(colorbar, fconfig)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_axis_off()
    ax.set_facecolor([0, 0, 0, 1])
    fig.tight_layout()
    return fig, ax
