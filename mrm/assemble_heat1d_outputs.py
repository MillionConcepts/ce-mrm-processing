"""
This script assembles the outputs of individual runs of heat1d.Model into a
table giving heatflow dependence on lat, ltst, h-parameter, albedo, depth,
thickness, and density. Must be run after execute_heat1d.py; generates inputs
for assemble_brightness_table.py. Various behaviors can be controlled by
modifying constants (described in code body below).
"""

from pathlib import Path
import pickle

import numpy as np
from numpy import dtype
import pandas as pd
import pyarrow as pa
from pyarrow import parquet as pq

from mrm.shared.console import print_inline

# data type specifications for columns of output parquet file
COLSPECS = {
    'temp': dtype('float32'),
    'z': dtype('float32'),
    'dz': dtype('float32'),
    'ltst': dtype('float32'),
    'rho': dtype('float32'),
    'albedo': dtype('float32'),
    'lat': dtype('float32'),
    'h': dtype('float32'),
    'layer_scale_factor': dtype('uint8'),
    'nskin': dtype('uint8'),
    'nskinbot': dtype('uint8'),
    'nyearseq': dtype('uint8'),
    'nperday': dtype('uint16'),
    'chi': dtype('float32'),
    'fourier_mesh_number': dtype('float32'),
    'dtsurf': dtype('float32'),
    'ztscale': dtype('float32'),
    'eccentricity': dtype('float32'),
    'obliquity': dtype('float32'),
    'qb': dtype('float32'),
    'ks': dtype('float32'),
    'kd': dtype('float32')
}

# where are the outputs of execute_heat1d.py, and where should we write the
# intermediate outputs of this script?
INPATH = Path(__file__).parent.parent / "pipeline_data/heat1d_chunks/"
OUTPATH = Path(__file__).parent.parent / "pipeline_data"
# Parameters we wish to vary in the output table. Must match parameters that
# varied across the heat1d run; will throw an error if it looks like they don't
# match.
VARIABLE_METADATA_FIELDS = ('qb', 'h', 'albedo', 'lat')
# time resolution of output table
TIME_SAMPLES = 240

if __name__ == "__main__":
    # TODO, maybe: something other than a bunch of pickle files would be better
    frames, invariants = [], {}
    files = tuple(filter(lambda p: p.suffix == '.pkl', INPATH.iterdir()))
    print(f"Processing {len(files)} model chunks...")
    for ix, file in enumerate(files):
        raw = pickle.load(file.open('rb'))
        temp = raw['data']['temp']
        layers = np.split(temp, temp.shape[0])
        time_rectified = [
            np.interp(
                np.linspace(0, temp.shape[1], TIME_SAMPLES),
                np.arange(temp.shape[1]),
                layer[0]
            )
            for layer in layers
        ]
        columns = {
            'temp': np.concatenate(time_rectified),
            # NOTE: model 0 time is noon; shift it to midnight
            'ltst': np.tile(
                (np.arange(0, 24, 24/TIME_SAMPLES) + 12) % 24, temp.shape[0]
            ),
            'z': np.repeat(raw['data']['z'], TIME_SAMPLES).ravel(),
            'dz': np.repeat(raw['data']['dz'], TIME_SAMPLES).ravel(),
            'rho': np.repeat(raw['data']['rho'], TIME_SAMPLES).ravel()
        }
        df = pd.DataFrame(columns)
        raw['meta']['lat'] = np.degrees(raw['meta']['lat'])
        for k, v in raw['meta'].items():
            if k in VARIABLE_METADATA_FIELDS:
                df[k] = v
            elif k not in invariants:
                invariants[k] = v
            elif invariants[k] != v:
                raise ValueError(
                    f"{k} was not listed in VARIABLE_METADATA_FIELDS, but its "
                    f"value varies across model outputs. This may indicate "
                    f"mismatched input products."
                )
        frames.append(
            df.astype({k: v for k, v in COLSPECS.items() if k in df.columns})
        )
        print_inline(f"processed {ix} chunks...")
    print("writing output table...")
    table = pa.Table.from_pandas(pd.concat(frames).reset_index(drop=True))
    del frames
    table = table.replace_schema_metadata(
        {k: str(v).encode('utf-8') for k, v in invariants.items()}
    )
    pq.write_table(
        table,
        OUTPATH / 'heat1d_table.parquet',
        row_group_size=500000,
        version="2.6"
    )
    print("Done.")
