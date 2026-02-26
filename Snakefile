# ============================================================
# RNA Flanking Sequence Pipeline — Main Snakefile
# ============================================================
# Stages:
#   1. parse_ids       → classify transcript IDs by source DB
#   2. resolve_ids     → map transcript IDs to gene + genomic coordinates
#   3. download_asm    → fetch & index genome assemblies (cached)
#   4. extract_seqs    → extract gene + flanking sequences
#   5. report          → generate summary HTML + TSV reports
# ============================================================

from pathlib import Path

configfile: "config/config.yaml"

# ── Convenience variables ────────────────────────────────────
RESULTS     = config["results_dir"]
LOGS        = config["logs_dir"]
BENCHMARKS  = config["benchmarks_dir"]
CACHE       = config["cache_dir"]
UPSTREAM    = config["upstream_bp"]
DOWNSTREAM  = config["downstream_bp"]

# ── Include rule modules ─────────────────────────────────────
include: "rules/parse_ids.smk"
include: "rules/resolve_ids.smk"
include: "rules/download_assemblies.smk"
include: "rules/extract_sequences.smk"
include: "rules/report.smk"

# ── Target rule ──────────────────────────────────────────────
rule all:
    input:
        fasta   = f"{RESULTS}/output.fasta",
        bed     = f"{RESULTS}/output.bed",
        report  = f"{RESULTS}/report.html",
        unresolv= f"{RESULTS}/unresolved.tsv",
        ambig   = f"{RESULTS}/ambiguous.tsv",
