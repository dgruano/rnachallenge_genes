"""
scripts/discover_config_db_urls.py
Stage 2 (config-backed DBs) — URL Table Discovery
===================================================

Reads a top-level config key (e.g., "plant_gtf_sources") and emits a
per-assembly URL table.  One row per species entry in the YAML.

Each row carries:
  db_source          — Snakemake params.db_source_override, or the YAML
                       entry key (e.g., "wormbase", "flybase")
  assembly_name      — human-readable assembly ID from config
  assembly_accession — GCF_/GCA_ if known in config; otherwise null
  fasta_url          — direct FASTA download URL if configured; else null
  gtf_url            — direct GTF/GFF download URL (the "url" config field)
  gtf_format         — "gtf" | "gff3" | null
  organism           — species/organism identifier from config

This output is consumed by merge_resolved (Stage 3) via a LEFT JOIN on
(assembly_name, db_source) to fill URL columns for config-backed DBs.

Params
------
snakemake.params.config_key : str
    Top-level key in snakemake.config, e.g. "plant_gtf_sources"
snakemake.params.db_source_override : str, optional
    If set, override the db_source for all rows (e.g., "ensembl_plants").
    If empty / not set, use the YAML entry key as db_source.

Output
------
snakemake.output.urls : path to TSV
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

URLS_COLS = [
    "db_source",
    "assembly_name",
    "assembly_accession",
    "fasta_url",
    "gtf_url",
    "gtf_format",
    "organism",
]

# ── Snakemake interface ──────────────────────────────────────────────────────
log = get_logger("discover_config_db_urls", snakemake.log[0])
out_urls = snakemake.output.urls
config_key = snakemake.params.config_key
db_source_override = getattr(snakemake.params, "db_source_override", "") or ""

# ── Main ─────────────────────────────────────────────────────────────────────
log.info(f"discover_config_db_urls: reading config key '{config_key}'")

sources = snakemake.config.get(config_key, {})
if not sources:
    log.warning(f"  Config key '{config_key}' is empty or missing; writing empty URL table")
    pd.DataFrame(columns=URLS_COLS).to_csv(out_urls, sep="\t", index=False)
    log.info("Done.")
    exit(0)

rows = []
for entry_key, entry in sources.items():
    if not isinstance(entry, dict):
        continue

    db_source = db_source_override if db_source_override else entry_key

    # 'url' is the GTF URL in all existing config files
    gtf_url = entry.get("url") or entry.get("gtf_url")
    assembly_accession = entry.get("assembly_accession")
    if assembly_accession is not None and str(assembly_accession).lower() in ("null", "none", ""):
        assembly_accession = pd.NA

    rows.append({
        "db_source": db_source,
        "assembly_name": entry.get("assembly_name", entry_key),
        "assembly_accession": assembly_accession,
        "fasta_url": entry.get("fasta_url"),
        "gtf_url": gtf_url,
        "gtf_format": entry.get("gtf_format"),
        "organism": entry.get("organism", entry_key),
    })

df_urls = pd.DataFrame(rows, columns=URLS_COLS)
log.info(f"Writing {len(df_urls)} assembly URL row(s) to {out_urls}")
df_urls.to_csv(out_urls, sep="\t", index=False)
log.info("Done.")
