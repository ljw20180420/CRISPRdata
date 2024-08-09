#!/usr/bin/env python

import argparse
import pathlib
import re
import numpy as np
import polars as pl
import tempfile
import os

# command line inputs
parser = argparse.ArgumentParser(description="preprocess alignments")
parser.add_argument("--data_dir", type=pathlib.Path, help="directory containing alignment files")
parser.add_argument("--score_quantile", type=float, default=0.05, help="alignment score quantile to filter score")
parser.add_argument("--min_score", type=float, default=-np.inf, help="min alignment score threshold")
parser.add_argument("--name", type=str, default="dataset.json", help="output name")
args = parser.parse_args()

# read data
from datasets import load_dataset, Features, Value
alg_features = Features({
    'index': Value('uint64'),
    'count': Value('uint64'),
    'score': Value('float32'),
    'ref_id': Value('uint32'),
    'up_dangle': Value('string'),
    'ref1_start': Value('int16'),
    'query1_start': Value('int16'),
    'ref1_end': Value('int16'),
    'query1_end': Value('int16'),
    'random_insert': Value('string'),
    'ref2_start': Value('int16'),
    'query2_start': Value('int16'),
    'ref2_end': Value('int16'),
    'query2_end': Value('int16'),
    'down_dangle': Value('string'),
    'cut1': Value('int16'),
    'cut2': Value('int16'),
    'ref': Value('string'),
    'query': Value('string')
})

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

# keep only count, score, ref1_end, ref2_start, cut1, cut2, ref1, ref2, random_insert
alg_ds = (
    load_dataset("csv", data_files=(args.data_dir / "*").as_posix(), num_proc=12, delimiter="\t", column_names=['index', 'count', 'score', 'ref_id', 'up_dangle', 'ref1_start', 'query1_start', 'ref1_end', 'query1_end', 'random_insert', 'ref2_start', 'query2_start', 'ref2_end', 'query2_end', 'down_dangle', 'cut1', 'cut2', 'ref', 'query'], features=alg_features, keep_default_na=False)
    # load_dataset("csv", data_files=(args.data_dir / "A-CtIPb-2.fq.alg.gz").as_posix(), num_proc=12, delimiter="\t", column_names=['index', 'count', 'score', 'ref_id', 'up_dangle', 'ref1_start', 'query1_start', 'ref1_end', 'query1_end', 'random_insert', 'ref2_start', 'query2_start', 'ref2_end', 'query2_end', 'down_dangle', 'cut1', 'cut2', 'ref', 'query'], features=alg_features, keep_default_na=False)
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
(
    pl.scan_parquet(temp_file, low_memory=True)
    .group_by(["ref1_end", "ref2_start", "random_insert", "cut1", "cut2", "ref1", "ref2"])
    .agg(pl.col("count").sum())
    .group_by(["cut1", "cut2", "ref1", "ref2"])
    .agg(pl.col("ref1_end"), pl.col("ref2_start"), pl.col("random_insert"), pl.col("count"))
    .group_by(["ref1", "ref2"])
    .agg(pl.col("cut1"), pl.col("cut2"), pl.col("ref1_end"), pl.col("ref2_start"), pl.col("random_insert"), pl.col("count"))
    .collect(streaming=True)
    .write_ndjson((args.data_dir.parent / args.name).as_posix())
    # .write_parquet((args.data_dir.parent / "dataset.parquet").as_posix())
)
os.remove(temp_file)