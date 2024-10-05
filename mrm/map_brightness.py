"""
This script produces the "tbmod" map files. Various aspects of its behavior
can be controlled by modifying constants (specifics described in comments in
the body of the code).

It requires the intermediate model files produced by running execute_heat1d.py,
assemble_heat1d_outputs.py, and assemble_brightness_table.py in order, along
with the following input files in pipeline_data/source_maps:

H-parameter-64ppd.fits
SLDEM2015_64_AZ_90S_90N_000_360.IMG
SLDEM2015_64_AZ_90S_90N_000_360.LBL
SLDEM2015_64_SL_90S_90N_000_360.IMG
SLDEM2015_64_SL_90S_90N_000_360.LBL
wac_hapke_604nm_70N70S_64ppd.fits
"""
from multiprocessing import Pool
from pathlib import Path
from typing import Literal

from mrm.processing.brightness_mapping import (
    make_filepaths, run_channel, prep_source_data, INT_PATHS, SRCMAP_FNS
)

# What are the root directories for source and intermediate data?
# Note that subdirectories are defined in make_filepaths().
PIPELINE_DATA_PATH = Path(__file__).parent.parent / "pipeline_data"
TBMOD_DATA_PATH = Path(__file__).parent.parent / "bundle/miscellaneous/tbmod/"

# which chanels should we run?
CHANNELS: tuple[Literal[1, 2, 3, 4], ...] = (1, 2, 3, 4)
# At what deg/pix did we generate the temp maps?
RESOLUTION: float = 1 / 32
# size of local time slices
TIMEBIN_SIZE: int = 2
# overwrite already-produced intermediate maps?
OVERWRITE_INTERMEDIATE_MAPS: bool = False
# maximum precision of derived LTST offsets (decimal places)
LTST_PRECISION = 1
# maximum absolute value of derived LTST offsets
MAX_TIME_CORRECTION: float = 1.5

# Fixed loss tangent parameter. Needs to be a value actually modeled when
# creating the table (in assemble_brightness_table.py).
C_CONSTANT = 2.5
# NOTE: drop this to 1 (or just remove threading, although there's no real
# overhead here) to run in 32 GB of RAM. Don't run this script at all in 16.
THREADS = 2

if __name__ == "__main__":
    # makes directories and also initializes globals in brightness_mapping.
    # should maybe be cleaner.
    make_filepaths(PIPELINE_DATA_PATH, TBMOD_DATA_PATH)
    timebins = [
        (start, start + TIMEBIN_SIZE)
        for start in range(0, 24, TIMEBIN_SIZE)
    ]
    if (
        OVERWRITE_INTERMEDIATE_MAPS is True
        or not all(INT_PATHS[k].exists() for k in SRCMAP_FNS.keys())
    ):
        print("(re)generating intermediate maps from source...")
        prep_source_data(RESOLUTION, True)
    pool = Pool(THREADS)
    kwargs = {
        'timebins': timebins,
        'resolution': RESOLUTION,
        'max_time_correction': MAX_TIME_CORRECTION,
        'overwrite_intermediate_maps': False,
        'ltst_precision': LTST_PRECISION,
        'c_constant': C_CONSTANT
    }
    for channel in CHANNELS:
        pool.apply_async(run_channel, kwds=(kwargs | {'channel': channel}))
    pool.close()
    pool.join()
