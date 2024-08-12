#!/usr/bin/env python

import argparse
import re
import numpy as np
import polars as pl
import tempfile
import os
import gzip
from datasets import load_dataset, Features, Value

# command line inputs
parser = argparse.ArgumentParser(description="preprocess alignments")
parser.add_argument("--score_quantile", type=float, default=0.05, help="alignment score quantile to filter score")
parser.add_argument("--min_score", type=float, default=-np.inf, help="min alignment score threshold")
args = parser.parse_args()

# ref1 and ref2 are splited by a lower case acgtn
acgtn_pattern = re.compile('[acgtn]')
def get_ref1_len(examples):
    refs = [ref.replace("-", "") for ref in examples['ref']]
    ref1_lens = [acgtn_pattern.search(ref[1:]).span()[1] + 1 for ref in refs]
    ref2_starts = [ref2_start - ref1_len for ref2_start, ref1_len in zip(examples['ref2_start'], ref1_lens)]
    cut2s = [cut2 - ref1_len for cut2, ref1_len in zip(examples['cut2'], ref1_lens)]
    ref1s = [ref[:ref1_len].upper() for ref1_len, ref in zip(ref1_lens, refs)]
    ref2s = [ref[ref1_len:].upper() for ref1_len, ref in zip(ref1_lens, refs)]
    return {
        "ref2_start": ref2_starts,
        "cut2": cut2s,
        "ref1": ref1s,
        "ref2": ref2s
    }

data_files = []
for author in os.listdir('data'):
    if not os.path.isdir(f'data/{author}/algs'):
        continue
    for data_file in os.listdir(f'data/{author}/algs'):
        data_files.append(f'data/{author}/algs/{data_file}')
# data_files = ["data/LE/algs/LAMSX-20230512-DCK-D-D-R1_L2_Q0421W0322.R2.fastq.alg.gz", "data/LE/algs/LAMSX-20230512-HPRT-D-D-R2_L2_Q0427W5344.R2.fastq.alg.gz"]

# keep only count, score, ref1_end, ref2_start, cut1, cut2, ref1, ref2, random_insert, author, file_name
alg_ds = (
    load_dataset(
        "csv",
        data_files=data_files,
        num_proc=12,
        delimiter="\t",
        column_names=['index', 'count', 'score', 'ref_id', 'up_dangle', 'ref1_start', 'query1_start', 'ref1_end', 'query1_end', 'random_insert', 'ref2_start', 'query2_start', 'ref2_end', 'query2_end', 'down_dangle', 'cut1', 'cut2', 'ref', 'query', 'author', 'file_name'], features=Features({
            'index': Value('uint64'),
            'count': Value('uint64'),
            'score': Value('float32'),
            'ref_id': Value('uint32'),
            'up_dangle': Value('string'),
            'ref1_start': Value('uint16'),
            'query1_start': Value('uint16'),
            'ref1_end': Value('uint16'),
            'query1_end': Value('uint16'),
            'random_insert': Value('string'),
            'ref2_start': Value('uint16'),
            'query2_start': Value('uint16'),
            'ref2_end': Value('uint16'),
            'query2_end': Value('uint16'),
            'down_dangle': Value('string'),
            'cut1': Value('uint16'),
            'cut2': Value('uint16'),
            'ref': Value('string'),
            'query': Value('string'),
            'author': Value('string'),
            'file_name': Value("string")
        }),
        keep_default_na=False
    )
    .remove_columns(['index', 'ref_id', 'up_dangle', 'ref1_start', 'query1_start', 'query1_end', 'query2_start', 'ref2_end', 'query2_end', 'down_dangle', 'query'])
    .map(get_ref1_len, batched=True)
    .remove_columns(["ref"])
)

# filter records by alignment score
score_thres = max(
    args.min_score,
    np.quantile(alg_ds['train']['score'], args.score_quantile)
)
alg_ds = (
    alg_ds.filter(lambda examples: [score >= score_thres for score in examples["score"]], batched=True)
    .remove_columns(["score"])['train']
)

# save results to a temperary parquet file
temp_file = tempfile.mkstemp(dir = "/home/ljw/sdc1/tmp")[1]
alg_ds.to_parquet(temp_file)
del alg_ds

# load parquet file into polars data frame and group by ...
# write the result to newline delimited JSON representation (parquet containing large_list item is not supported by huggingface datasets at this time)
with gzip.open("dataset.json.gz", 'wb') as gfd:
    (
        pl.scan_parquet(temp_file, low_memory=True)
        .group_by(["ref1_end", "ref2_start", "random_insert", "cut1", "cut2", "ref1", "ref2", "file_name", "author"])
        .agg(pl.col("count").sum())
        .group_by(["cut1", "cut2", "ref1", "ref2", "file_name", "author"])
        .agg(pl.col("ref1_end"), pl.col("ref2_start"), pl.col("random_insert"), pl.col("count"))
        .group_by(["ref1", "ref2", "file_name", "author"])
        .agg(pl.col("cut1"), pl.col("cut2"), pl.col("ref1_end"), pl.col("ref2_start"), pl.col("random_insert"), pl.col("count"))
        .collect(streaming=True)
        .write_ndjson(gfd)
    )
os.remove(temp_file)