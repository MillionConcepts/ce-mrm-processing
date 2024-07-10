"""to be run after assemble_heat1d_outputs"""

import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pyarrow as pa
from pyarrow import parquet as pq

from mrm.shared.constants import CHANNEL_FREQ
from mrm.modeling.brightness import simulate_brightnesses
from mrm.shared.pqutils import pq_metadict

INPUT_TABLE = (
    Path(__file__).parent.parent / "pipeline_data/heat1d_table.parquet"
)
OUTPUT_TABLE = Path(__file__).parent.parent / "pipeline_data/tb_table.parquet"
FREQUENCIES = tuple(CHANNEL_FREQ.values())
A, B, TIO2_PERCENTAGE = 0, 0.312, 0
C_VALUES = np.arange(1.5, 3.51, 0.02, dtype='f4').round(3)
# NOTE: you'll want to drop THREADS to 5-6 if you're running with 32 GB of RAM,
#  and to 2-3 if running in 16 (and 16 may be a bit sketchy at all).
THREADS = 12

if __name__ == "__main__":
    meta = pq_metadict(INPUT_TABLE) | {
        'a': str(A).encode('utf-8'),
        'b': str(B).encode('utf-8'),
        'tio2_percentage': str(TIO2_PERCENTAGE).encode('utf-8')
    }
    df = pq.read_table(INPUT_TABLE).to_pandas()
    pool, table, results = Pool(THREADS), None, {}
    params = {
        # NOTE: this assumes LTST values are shared among all groups
        'ltst': sorted(df['ltst'].unique()),
        'frequencies': FREQUENCIES,
        'c_values': C_VALUES,
        'a': A,
        'b': B,
        'tio2_percentage': TIO2_PERCENTAGE
    }
    lats = list(df['lat'].unique())
    while (len(lats) > 0) or (len(results) > 0):
        if (len(lats) > 0) and (len(results) < THREADS):
            lat = lats.pop()
            results[lat] = pool.apply_async(
                simulate_brightnesses, (lat, df.loc[df['lat'] == lat]), params
            )
        elif len(lats) == 0:
            pool.close()
        completed = [k for k in results.keys() if results[k].ready()]
        for c in completed:
            if table is None:
                table = results.pop(c).get()
                print("created table")
            else:
                # NOTE: this clunky operation is a time-for-space trade
                table = pa.concat_tables(
                    [table, results.pop(c).get()]
                ).unify_dictionaries().combine_chunks()
                print(f"merged table ({len(lats) + len(results)} remaining)")
        time.sleep(0.25)
    pool.terminate()
    print("writing table")
    pq.write_table(
        table, OUTPUT_TABLE, row_group_size=1000000, version="2.6"
    )
    print("done")
