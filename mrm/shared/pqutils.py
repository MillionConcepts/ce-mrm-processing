from pyarrow import parquet


def pq_metadict(tablepath):
    parquet_meta = parquet.read_metadata(tablepath).metadata
    return {
        k: v for k, v in parquet_meta.items() if k != b'ARROW:schema'
    }
