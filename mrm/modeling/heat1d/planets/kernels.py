"""This module helps to manage SPICE kernels, incl. downloading and listing
loaded kernels (which unbelievably is not available from SPICE directly).
"""
from pathlib import Path

import spiceypy as spice
import requests

KERNELROOT = Path(__file__).parent / "kernels"
KERNELROOT.mkdir(exist_ok=True)

download_root = Path("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/")

generic_kernel_names = [
    "lsk/naif0012.tls",
    "pck/pck00011.tpc",
    "pck/de-403-masses.tpc",
    "spk/planets/de440.bsp",
    "spk/satellites/mar097.bsp",
]

generic_kernels_paths = [KERNELROOT / i for i in generic_kernel_names]


def do_download(source, target):
    """Download source url to target path.

    Parameters
    ----------
    source: url <str>
    target: pathlib.Path
    """
    target.parent.mkdir(parents=True, exist_ok=True)
    response = requests.get(source)
    response.raise_for_status()
    with open(target, "wb") as stream:
        stream.write(response.content)


def download_generic_kernels(kernel=None):
    "Download all kernels as required by generic_kernel_list."
    if kernel is None:
        dl_urls = [download_root / i for i in generic_kernel_names]
        for dl_url, savepath in zip(dl_urls, generic_kernels_paths):
            do_download(dl_url, savepath)
    else:
        dl_url = download_root / kernel.relative_to(KERNELROOT)
        print(f"Downloading\n{dl_url}\ninto\n{kernel}\n")
        do_download(dl_url, kernel)


def check_generic_kernels():
    "Check for existence of generic_kernels and download if not there."
    for kernel in generic_kernels_paths:
        if not kernel.exists():
            print(f"Cannot find generic kernel {kernel.name} locally.")
            download_generic_kernels(kernel)


def load_generic_kernels():
    """Load all kernels in generic_kernels list.

    Loads pure planetary bodies meta-kernel without spacecraft data.

    Calls `check_generic_kernels()` which downloads the kernels if they are
    not there.
    """
    check_generic_kernels()
    for kernel in generic_kernels_paths:
        spice.furnsh(str(kernel))


def load_planet_masses_kernel():
    spice.furnsh(KERNELROOT / "pck/de-403-masses.tpc")


def show_loaded_kernels():
    "Print overview of loaded kernels."
    count = spice.ktotal("all")
    if count == 0:
        print("No kernels loaded at this time.")
    else:
        print("The loaded files are:\n(paths relative to kernels.KERNELROOT)\n")
    for which in range(count):
        print(which)
        out = spice.kdata(which, "all", 100, 100, 100)
        print("Position:", which)
        p = Path(out[0])
        print("Path", p.relative_to(KERNELROOT))
        print("Type:", out[1])
        print("Source:", out[2])
        print("Handle:", out[3])
        print("Found:", out[4])
