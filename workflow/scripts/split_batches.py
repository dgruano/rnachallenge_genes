"""
scripts/split_batches.py

Split resolved_ids.tsv into N batch TSVs and write batch_manifest.txt.

Snakemake interface:
    snakemake.input.resolved      — results/resolved_ids.tsv
    snakemake.output.manifest     — results/batch_manifest.txt
    snakemake.output.batch_dir    — results/batches/ (directory)
    snakemake.params.batch_size
    snakemake.log[0]
"""

import sys
from math import ceil
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

log = get_logger("split_batches", snakemake.log[0])

resolved_path = snakemake.input.resolved
manifest_path = snakemake.output.manifest
batch_dir = Path(snakemake.output.batch_dir)
batch_size = int(snakemake.params.batch_size)

df = pd.read_csv(resolved_path, sep="\t")
n_rows = len(df)
n_batches = ceil(n_rows / batch_size)

log.info(f"Splitting {n_rows} rows into {n_batches} batches of ~{batch_size}")

batch_dir.mkdir(parents=True, exist_ok=True)

batch_ids = []
for i in range(n_batches):
    batch_id = f"batch_{i:04d}"
    chunk = df.iloc[i * batch_size : (i + 1) * batch_size]
    out_path = batch_dir / f"{batch_id}.tsv"
    chunk.to_csv(out_path, sep="\t", index=False)
    batch_ids.append(batch_id)
    log.info(f"  Written {out_path} ({len(chunk)} rows)")

with open(manifest_path, "w") as fh:
    for batch_id in batch_ids:
        fh.write(f"{batch_id}\n")

log.info(f"Manifest written: {manifest_path} ({len(batch_ids)} batches)")
