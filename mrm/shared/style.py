from pathlib import Path

import matplotlib as mpl
from matplotlib import colors as mcolor
from matplotlib.cm import register_cmap
import matplotlib.pyplot as plt
import numpy as np
"""
style/plot control functions, along with colormaps. intended for browse 
products and optional figures.
"""
# NOTE: partially vendored from moonbow.


EASYDARK = mpl.RcParams(
    {
        "axes.facecolor": "black",
        "axes.labelcolor": (0, 0, 0, 0),
        "figure.edgecolor": (0, 0, 0, 0),
        "figure.facecolor": "black",
        "patch.edgecolor": (0, 0, 0, 0),
        "savefig.edgecolor": "black",
        "savefig.facecolor": "black",
        "text.color": "white",
        "xtick.color": "white",
        "ytick.color": "white",
    }
)
FONT_PATH = Path(__file__).parent.parent.parent / "static/fonts"


def make_orange_teal_cmap():
    teal = (98, 252, 232)
    orange = (255, 151, 41)
    half_len = 256
    vals = np.ones((half_len * 2, 4))
    vals[0:half_len, 0] = np.linspace(orange[0] / half_len, 0, half_len)
    vals[0:half_len, 1] = np.linspace(orange[1] / half_len, 0, half_len)
    vals[0:half_len, 2] = np.linspace(orange[2] / half_len, 0, half_len)
    vals[half_len:, 0] = np.linspace(0, teal[0] / half_len, half_len)
    vals[half_len:, 1] = np.linspace(0, teal[1] / half_len, half_len)
    vals[half_len:, 2] = np.linspace(0, teal[2] / half_len, half_len)
    return mcolor.ListedColormap(vals, name="orange_teal")


def make_sunset():
    return mcolor.LinearSegmentedColormap.from_list(
        'vaportrail',
        [
            "#000000",
            "#210740",
            "#590080",
            "#BF5AAD",
            "#FF96B6",
        ]
)


def make_vaportrail():
    return mcolor.LinearSegmentedColormap.from_list(
        'sunset',
        [
            "#000000",
            "#001F2E",
            "#2B3BB8",
            "#CD6BF5",
            "#FFAD96",
        ]
    )


def make_neon():
    return mcolor.LinearSegmentedColormap.from_list(
        'neon',
        [
            "#FF00D4",
            "#941212",
            "#000000",
            "#3B73F5",
            "#9CFCFF",
        ]
    )


def make_bowmoon():
    return mcolor.LinearSegmentedColormap.from_list(
        'bowmoon',
        [
            "#351D21",
            "#6F3812",
            "#766C07",
            "#10A737",
            "#96BEF6",
        ]
    )


def make_moonbow():
    return mcolor.LinearSegmentedColormap.from_list(
        'moonbow',
        [
            "#3E055C",
            "#462BBC",
            "#2866C0",
            "#28A72B",
            "#CFB675",
            "#F9DCFC",
        ]
    )


def register_cmaps():
    for func in (
        make_orange_teal_cmap,
        make_sunset,
        make_vaportrail,
        make_neon,
        make_moonbow,
        make_bowmoon
    ):
        try:
            cmap = func()
            register_cmap(name=cmap.name, cmap=cmap)
            register_cmap(name=f"{cmap.name}_r", cmap=cmap.reversed())
        except ValueError:  # most likely someone registered it already
            continue


class NoShow:
    """
    context manager to make matplotlib not do annoying interactive
    things and auto-close figures. if you want to get figures out,
    pass close=False
    """

    def __init__(self, close: bool = True):
        self.close = close

    def __enter__(self):
        plt.ioff()

    def __exit__(self, *args, **kwargs):
        if self.close is True:
            plt.close('all')
        plt.ion()
