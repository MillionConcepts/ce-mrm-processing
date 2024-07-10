"""
Handler script for physical temperature / heatflow modeling step. Generates
.pkl files, by default in pipeline_data/heat1d_chunks, that can be used as
inputs to assemble_heat1d_outputs.py. Various behaviors can be controlled by
modifying constants in the script.
"""

from pathlib import Path

import numpy as np

from mrm.modeling.heat1d.runners import run_parallel

PARAMS = {
    "eccentricity": 0,
    "obliquity": 0,
    # (0.017 unlabeled heatflow, probably mare?)
    # original C code uses 0.018.
    # "qb": 0.008,  # highlands
    # Solid (phonon) conductivity at surface [W.m-1.K-1]
    "qb": 0.013,  # compromise
    "ks": 8.0e-4,
    # Solid (phonon) conductivity at depth [W.m-1.K-1]
    "kd": 3.8e-3,
    "layer_scale_factor": 10,  # NN / N / n
    "nskinbot": 60,  # NSKINBOT / b
    "nskin": 20,  # NSKIN / m
    "nyearseq": 80,
    "nperday": 480,
}
"""
Constant heatflow model parameters. Used as arguments to every call to 
heat1d.Model.__init__().
"""

AXBOUNDS = {
    "latitudes": {
        "range": (0, 75),
        "interval": 2,
        "precision": 2,
    },
    "h_parameters": {
        "range": (0.02, 0.13),
        "interval": 0.01,
        "precision": 3,
    },
    "albedos": {
        "range": (0.02, 0.16),
        "interval": 0.01,
        "precision": 3,
    },
}
"""
Variable heatflow model parameters. This script instantiates and executes one
copy of heat1d.Model for each element of the cartesian product of the 
evenly-spaced arrays of values specified by these dicts.
"""

AXES = {
    ax: np.round(
        np.arange(
            AXBOUNDS[ax]["range"][0],
            AXBOUNDS[ax]["range"][1] + AXBOUNDS[ax]["interval"],
            AXBOUNDS[ax]["interval"],
        ),
        AXBOUNDS[ax]["precision"],
    )
    for ax in ("latitudes", "h_parameters", "albedos")
}
"""Parameter arrays assembled from AXBOUNDS."""

OUTPATH = Path(__file__).parent.parent / "pipeline_data" / "heat1d_chunks"
"""Where should we write model outputs?"""
SKIP_EXISTING = False
"""
If True, skip generating any permutation of the model for which an output file
already exists. Can be used to restart a crashed or paused run from where it 
left off.of it
"""
THREADS = 14
"""
How many copies of heat1d.Model should we run in parallel? Parallelism here is
constrained almost entirely by processor count; working memory is very 
unlikely to be a limiting factor at default settings (~150 MB per thread).  
"""

if __name__ == "__main__":
    OUTPATH.mkdir(exist_ok=True)
    run_parallel(
        PARAMS,
        **AXES,
        outpath=OUTPATH,
        skip_existing=SKIP_EXISTING,
        threads=THREADS
    )
