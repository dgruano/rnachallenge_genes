"""
scripts/aggregate_sequences.py

Concatenate per-batch FASTA, BED, and failed TSV files into final outputs.

FASTA: binary concatenation (FASTA files concatenate safely).
BED: header from first non-empty batch; append data rows from all others.
Failed TSV: same header-preservation logic; header-only file if all empty.

Snakemake interface:
    snakemake.input.fastas   — list of batch .fasta files
    snakemake.input.beds     — list of batch .bed files
    snakemake.input.faileds  — list of batch .failed.tsv files
    snakemake.output.fasta   — results/output.fasta
    snakemake.output.bed     — results/output.bed
    snakemake.output.failed  — results/extraction_failed.tsv
    snakemake.log[0]
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

log = get_logger("aggregate_sequences", snakemake.log[0])

fastas = snakemake.input.fastas
beds = snakemake.input.beds
faileds = snakemake.input.faileds
out_fasta = Path(snakemake.output.fasta)
out_bed = Path(snakemake.output.bed)
out_failed = Path(snakemake.output.failed)

FAILED_COLUMNS = [
    "transcript_id",
    "assembly_accession",
    "chrom",
    "db_source",
    "fail_reason",
]

# ── FASTA: binary concatenation ──────────────────────────────
log.info(f"Concatenating {len(fastas)} FASTA files → {out_fasta}")
out_fasta.parent.mkdir(parents=True, exist_ok=True)
with open(out_fasta, "wb") as fout:
    for path in fastas:
        fout.write(Path(path).read_bytes())

# ── BED: header from first non-empty file, then data rows ────
log.info(f"Concatenating {len(beds)} BED files → {out_bed}")
header_written = False
with open(out_bed, "w") as fout:
    for path in beds:
        lines = Path(path).read_text().splitlines(keepends=True)
        if not lines:
            continue
        if not header_written:
            fout.writelines(lines)  # header + data rows
            header_written = True
        else:
            fout.writelines(lines[1:])  # skip header, append data rows
if not header_written:
    out_bed.write_text("")

# ── Failed TSV: same header-preservation logic ───────────────
log.info(f"Concatenating {len(faileds)} failed TSV files → {out_failed}")
header_written = False
with open(out_failed, "w") as fout:
    for path in faileds:
        lines = Path(path).read_text().splitlines(keepends=True)
        if not lines:
            continue
        if not header_written:
            fout.writelines(lines)
            header_written = True
        else:
            fout.writelines(lines[1:])
if not header_written:
    with open(out_failed, "w") as fout:
        fout.write("\t".join(FAILED_COLUMNS) + "\n")

log.info("aggregate_sequences complete")
