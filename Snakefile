# ============================================================
# RNA Flanking Sequence Pipeline — Main Snakefile
# ============================================================
# DAG overview:
#
#   parse_ids
#     ├── resolve_ids (NCBI + UCSC)                       ─┐
#     └── detect_ensembl_species [checkpoint]              │
#               └── biomart_lookup (×N species) [wrapper]  ├── merge_resolved
#                         └── join_ensembl_results         ─┘
#                                                               │
#                                                     download_assemblies [checkpoint]
#                                                               │
#                                                     extract_sequences
#                                                               │
#                                                           report
# ============================================================

from pathlib import Path

configfile: "config/config.yaml"

RESULTS    = config["results_dir"]
LOGS       = config["logs_dir"]
BENCHMARKS = config["benchmarks_dir"]
CACHE      = config["cache_dir"]
UPSTREAM   = config["upstream_bp"]
DOWNSTREAM = config["downstream_bp"]

# ── Rule modules ─────────────────────────────────────────────
include: "rules/parse_ids.smk"
include: "rules/resolve_ids.smk"                # NCBI + UCSC
include: "rules/detect_ensembl_species.smk"     # checkpoint: infer species
include: "rules/biomart_lookup.smk"             # wrapper: per-species BioMart
include: "rules/join_ensembl_results.smk"       # join BioMart tables → Ensembl resolved
include: "rules/merge_resolved.smk"             # unify all three DB streams
include: "rules/download_assemblies.smk"        # checkpoint: cache genome FASTAs
include: "rules/extract_sequences.smk"
include: "rules/report.smk"

# ── Target rule ──────────────────────────────────────────────
rule all:
    input:
        fasta    = f"{RESULTS}/output.fasta",
        bed      = f"{RESULTS}/output.bed",
        report   = f"{RESULTS}/report.html",
        unresolv = f"{RESULTS}/unresolved.tsv",
        ambig    = f"{RESULTS}/ambiguous.tsv",
