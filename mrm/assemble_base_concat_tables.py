"""
This script concatenates the contents of the source data collections into
intermediate parquet files used by assemble_data_tables.py. It should be run
before running anything else. Input and output paths can be modified by editing
constants described in the body of the script.
"""

from pathlib import Path
import re

from cytoolz import first
from hostess.directory import index_breadth_first
from lhorizon.lhorizon_utils import utc_to_jd, utc_to_tdb
import pandas as pd
import pyarrow as pa
from pyarrow import parquet as pq

from mrm.shared.console import print_inline

ORBIT_PATTERN = re.compile(r".*_(\d{4})_[AB]\.2[BC]")
COLUMNS = [
    'time', 't1', 't2', 't3', 't4', 'theta', 'phi', 'lon', 'lat', 'd', 'flag'
]

# where are the L2C files stored?
L2C_PATH = Path(__file__).parent.parent / "bundle" / "data_source"
# where shall we write the intermediate table?
INT_PATH = Path(__file__).parent.parent / "pipeline_data" / "mrm_l2c_tables"
# which orbiter(s) shall we build intermediate table(s) for?
ORBITERS = ("ce1", "ce2")


def read_mrm_table(path):
    with open(path) as stream:
        for line in stream:
            if line.startswith("END\n"):
                break
        next(stream)
        table = pd.read_csv(stream, names=COLUMNS, sep=r'\s+')
    return table


def load_mrm_source_table(path):
    files = index_breadth_first(path)
    tables = {}
    for rec in files:
        if not (p := rec['path']).endswith("2C"):
            continue
        print_inline(f"loading {Path(p).name}...")
        tables[ORBIT_PATTERN.search(p).group(1)] = read_mrm_table(p)
    empties = {k for k, v in tables.items() if len(v) == 0}
    output = []
    print("organizing output...")
    while len(tables) > 0:
        orbit = first(tables.keys())
        table = tables.pop(orbit)
        if len(table) > 0:
            table['orbit'] = int(orbit)
            output.append(table)
    print("concatenating tables...")
    output = pd.concat(output).sort_values(by="orbit").reset_index(drop=True)
    return output, empties


def assemble_intermediate_table(orb):
    print("loading files...")
    l2c, _ = load_mrm_source_table(L2C_PATH / orb)
    l2c['time'] = l2c['time'].str.strip('Z').astype('datetime64[ms]')
    print("converting times...")
    l2c['jd'] = utc_to_jd(l2c['time'])
    l2c['tdb'] = utc_to_tdb(l2c['time'])
    l2c = l2c.sort_values(by=["orbit", "jd"]).reset_index(drop=True)
    arrays = {
        'time': l2c['time'].values,
        'jd': l2c['jd'].values,
        'tdb': l2c['tdb'].values.astype('datetime64[ms]'),
        't1': l2c['t1'].values.astype('float32'),
        't2': l2c['t2'].values.astype('float32'),
        't3': l2c['t3'].values.astype('float32'),
        't4': l2c['t4'].values.astype('float32'),
        'theta': l2c['theta'].values.astype('float32'),
        'phi': l2c['phi'].values.astype('float32'),
        'lat': l2c['lat'].values.astype('float32'),
        'lon': l2c['lon'].values.astype('float32'),
        'd': l2c['d'].values.astype('float32'),
        'flag': l2c['flag'].values.astype('str'),
        'orbit': l2c['orbit'].values.astype('int16')
    }
    del l2c
    l2c = pa.Table.from_arrays(list(arrays.values()), list(arrays.keys()))
    del arrays
    INT_PATH.mkdir(exist_ok=True, parents=True)
    print("writing intermediate parquet file...")
    pq.write_table(
        l2c, INT_PATH / f'{orb}_l2c_complete.parquet', version="2.6"
    )


if __name__ == "__main__":
    for o in ORBITERS:
        print(f"assembling table for {o}...")
        assemble_intermediate_table(o)
    print("done.")
