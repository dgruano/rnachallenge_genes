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
    # WormBase (full header-encoded IDs or gene-style IDs)
    (
        "wormbase",
        re.compile(
            r".*_wormbase:known_chromosome:WBcel\d+:",
            re.IGNORECASE,
        ),
    ),
    (
        "wormbase",
        re.compile(
            r"^[A-Z0-9]{1,3}\d+[A-Z]?\d*\.\d+[a-z]?(?:\.\d+)?$",
            re.IGNORECASE,
        ),
    ),
    # Plant gene/transcript IDs (Solanum, Oryza, Glycine, etc.)
    (
        "plant",
        re.compile(r"^Solyc\d+g\d+\.\d+\.\d+$", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^OS\d+T\d+(?:_\d+)?(?:_cdna)?$", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^Glyma\.\d{2}G\d{6}(?:\.\d+)?$", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^AT\dG\d+$", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^Zm\d+g\d+$", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^GRMZM\w+", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^LOC_Os\d+g\d+", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^Os\d+g\d+", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^TraesCS\w+\.\d+$", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^Bradi\dg\d+$", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^Bra\d+", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^BnaA\d+g\d+$", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^BnaC\d+g\d+$", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^Bo\dg\d+$", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^AET\w+", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^Amtr_\d+", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^evm\.model\.\w+\.\d+", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^Cre\d+\.g\d+", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^Pp\d+s\d+_\d+", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^Medtr\dg\d+", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^GSMUA_\w+", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^OB\d+g\d+", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^Si\d+g\d+", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^Thecc1EG\d+", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^orange1\.1g\d+", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^cassava\d+\.\d+", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^AC\d+\.\d+", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^Potri\.\d+G\d+$", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^Sobic\.\d+G\d+$", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^VIT_\d+s\d+$", re.IGNORECASE),
    ),
    (
        "plant",
        re.compile(r"^PGSC\d+DM\w+", re.IGNORECASE),
    ),
    # FlyBase (Drosophila)
    (
        "flybase",
        re.compile(r"^FBtr\d+$", re.IGNORECASE),
    ),
    (
        "flybase",
        re.compile(r"^FBgn\d+$", re.IGNORECASE),
    ),
    # WormBase gene IDs
    (
        "wormbase",
        re.compile(r"^WBGene\d+$", re.IGNORECASE),
    ),
    # SGD (yeast)
    (
        "sgd",
        re.compile(r"^Y[A-P][LR]\d+[WC](?:_[A-Z])?$", re.IGNORECASE),
    ),
    (
        "sgd",
        re.compile(r"^Q\d{4}$", re.IGNORECASE),
    ),
    (
        "sgd",
        re.compile(r"^Source:SGD;Acc:S\d+$", re.IGNORECASE),
    ),
]

# Embedded accessions found inside longer headers (GI/RefSeq, etc.)
EMBEDDED_PATTERNS: list[tuple[str, re.Pattern]] = [
    (
        "ncbi",
        re.compile(
            r"(?:ref[\|_])?((?:NM|NR|XM|XR|NP|XP|NG|NC|NT|NW|NZ)_\d+(?:\.\d+)?)",
            re.IGNORECASE,
        ),
    ),
    (
        "ensembl",
        re.compile(r"(ENS[A-Z]*T\d{11}(?:\.\d+)?)", re.IGNORECASE),
    ),
    (
        "ucsc",
        re.compile(r"(uc\d{3}[a-z]{3}\.\d+)", re.IGNORECASE),
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


def find_embedded_accession(header: str) -> tuple[str, str] | None:
    """Return (accession, db) if a known ID is embedded in the header."""
    for db, pattern in EMBEDDED_PATTERNS:
        match = pattern.search(header)
        if match:
            accession = match.group(1)
            return accession, db
    return None


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

        if db_source == "unknown":
            embedded = find_embedded_accession(raw_header)
            if embedded:
                transcript_id, db_source = embedded
                log.debug(
                    f"  Fallback extraction: {transcript_id!r} → {db_source} "
                    f"(from header: {raw_header!r})"
                )

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
