"""
scripts/parse_ids.py
Stage 1 — Parse & Classify Transcript IDs
==========================================
Reads one or more input FASTA files (via snakemake.input.fastas),
extracts the transcript ID from each header, and classifies each
ID into one of: ncbi | ensembl | ucsc | unknown.

Outputs
-------
classified_ids.tsv  : transcript_id, db_source, raw_header, source_file
unresolved.tsv      : transcript_id, raw_header, source_file, reason
"""

import re
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("parse_ids", snakemake.log[0])
fastas = snakemake.input.fastas
out_cls = snakemake.output.classified
out_unk = snakemake.output.unknown

# ── ID classification patterns ───────────────────────────────
# Order matters: more specific patterns first.
DB_PATTERNS: list[tuple[str, re.Pattern]] = [
    # NCBI RefSeq mRNA / ncRNA / predicted
    (
        "ncbi",
        re.compile(
            r"^(NM|NR|XM|XR|NP|XP|NG|NC|NT|NW|NZ)_\d+(\.\d+)?$",
            re.IGNORECASE,
        ),
    ),
    # Ensembl — covers human (ENST), mouse (ENSMUST), rat (ENSRNOT), etc.
    (
        "ensembl",
        re.compile(
            r"^ENS[A-Z]*T\d{11}(\.\d+)?$",
            re.IGNORECASE,
        ),
    ),
    # UCSC — e.g. uc001aaa.3
    (
        "ucsc",
        re.compile(
            r"^uc[0-9]{3}[a-z]{3}\.\d+$",
            re.IGNORECASE,
        ),
    ),
]


def classify_id(transcript_id: str) -> str:
    """Return the database name for a transcript ID, or 'unknown'."""
    for db, pattern in DB_PATTERNS:
        if pattern.match(transcript_id):
            return db
    return "unknown"


def extract_transcript_id(header: str) -> str:
    """
    Extract the primary accession from a FASTA header.
    Handles formats like:
      >NM_001234.3 Homo sapiens ...
      >ENST00000123456.7 ...
      >uc001aaa.3 ...
      >NM_001234.3|gene=BRCA1|...
    """
    # Strip leading '>'
    header = header.lstrip(">").strip()
    # Take first whitespace-delimited token, then first pipe-delimited part
    token = header.split()[0].split("|")[0]
    return token


# ── Main ─────────────────────────────────────────────────────
log.info("Stage 1: Parsing and classifying transcript IDs")
log.info(f"Input FASTA files: {fastas}")

classified_rows: list[dict] = []
unknown_rows: list[dict] = []

total_records = 0

for fasta_path in fastas:
    fasta_path = str(fasta_path)
    log.info(f"Processing: {fasta_path}")
    try:
        records = list(SeqIO.parse(fasta_path, "fasta"))
    except Exception as exc:
        log.error(f"Failed to parse {fasta_path}: {exc}")
        continue

    log.info(f"  Found {len(records)} records in {fasta_path}")
    total_records += len(records)

    for rec in records:
        raw_header = rec.description
        transcript_id = extract_transcript_id(rec.id)
        db_source = classify_id(transcript_id)

        row = {
            "transcript_id": transcript_id,
            "raw_header": raw_header,
            "source_file": fasta_path,
        }

        if db_source == "unknown":
            log.debug(
                f"  Unknown ID format: {transcript_id!r} (header: {raw_header!r})"
            )
            unknown_rows.append(
                {**row, "reason": "ID format not recognised as NCBI/Ensembl/UCSC"}
            )
        else:
            log.debug(f"  Classified {transcript_id!r} → {db_source}")
            classified_rows.append({**row, "db_source": db_source})

# ── Write outputs ─────────────────────────────────────────────
df_cls = pd.DataFrame(
    classified_rows, columns=["transcript_id", "db_source", "raw_header", "source_file"]
)
df_unk = pd.DataFrame(
    unknown_rows, columns=["transcript_id", "raw_header", "source_file", "reason"]
)

# De-duplicate by transcript_id (same ID in multiple files → keep first occurrence)
before = len(df_cls)
df_cls = df_cls.drop_duplicates(subset="transcript_id", keep="first")
if before != len(df_cls):
    log.warning(
        f"Removed {before - len(df_cls)} duplicate transcript IDs (kept first occurrence)"
    )

df_cls.to_csv(out_cls, sep="\t", index=False)
df_unk.to_csv(out_unk, sep="\t", index=False)

# ── Summary log ──────────────────────────────────────────────
log.info("=" * 60)
log.info(f"Total records parsed         : {total_records}")
log.info(f"Classified (NCBI)            : {(df_cls['db_source'] == 'ncbi').sum()}")
log.info(f"Classified (Ensembl)         : {(df_cls['db_source'] == 'ensembl').sum()}")
log.info(f"Classified (UCSC)            : {(df_cls['db_source'] == 'ucsc').sum()}")
log.info(f"Unknown / unclassified       : {len(df_unk)}")
log.info(f"Written classified  → {out_cls}")
log.info(f"Written unresolved  → {out_unk}")
log.info("Stage 1 complete.")
