# ============================================================
# RNA Flanking Sequence Pipeline — Main Snakefile
# ============================================================
# DAG overview:
#
#   parse_ids
#     ├── resolve_ids (NCBI + UCSC)                       ─┐
#     └── detect_ensembl_species [checkpoint]              │
#               └── biomart_lookup (×N species) [wrapper]  ├── merge_resolved
#                         └── join_ensembl_results        ─┘
#                                                               │
#                                                     download_assemblies [checkpoint]
#                                                               │
#                                                     extract_sequences
#                                                               │
#                                                           report
# ============================================================

from pathlib import Path

configfile: "config/config.yaml"
configfile: "config/plant_gtf_sources.yaml"
configfile: "config/tool_sources.yaml"

RESULTS    = config["results_dir"]
LOGS       = config["logs_dir"]
BENCHMARKS = config["benchmarks_dir"]
CACHE      = config["cache_dir"]
UPSTREAM   = config["upstream_bp"]
DOWNSTREAM = config["downstream_bp"]

# ── Rule modules ─────────────────────────────────────────────
include: "workflow/rules/source_mapping.smk"          # Stage 0: download tool FASTAs + map IDs
include: "workflow/rules/parse_ids.smk"
include: "workflow/rules/download_metadata.smk"          # metadata table downloads
include: "workflow/rules/resolve_ids.smk"                # NCBI + UCSC
include: "workflow/rules/resolve_external_ids.smk"       # plants + WormBase
include: "workflow/rules/biomart_plant_batch.smk"        # BioMart batch for plants
include: "workflow/rules/resolve_plant_gtf.smk"          # GTF fallback for BioMart failures
include: "workflow/rules/gramene_resolver.smk"           # Gramene API for legacy IDs
include: "workflow/rules/detect_ensembl_species.smk"     # checkpoint: infer species
include: "workflow/rules/biomart_lookup.smk"             # wrapper: per-species BioMart
include: "workflow/rules/join_ensembl_results.smk"       # join BioMart tables → Ensembl resolved
include: "workflow/rules/resolve_noncode.smk"            # NONCODE v5 transcript/gene IDs
include: "workflow/rules/resolve_noncode_v4.smk"         # NONCODEv4 fallback for v5-unresolved IDs
include: "workflow/rules/resolve_noncode_2016.smk"       # NONCODE2016 fallback (existence only)
include: "workflow/rules/resolve_ncbi_genbank.smk"       # EPost→EFetch second-pass resolver
include: "workflow/rules/resolve_abandoned_accessions.smk"  # GTF-based third-pass resolver
include: "workflow/rules/merge_resolved.smk"             # unify all three DB streams
include: "workflow/rules/download_assemblies.smk"        # checkpoint: cache genome FASTAs
include: "workflow/rules/extract_sequences.smk"
include: "workflow/rules/report.smk"

# ── Target rule ──────────────────────────────────────────────
rule all:
    input:
        # Stage 0 — tool source mapping (preprocessing)
        tool_map      = f"{RESULTS}/tool_source_map.tsv",
        tool_stats    = f"{RESULTS}/tool_source_stats.tsv",
        tool_unmatched = f"{RESULTS}/tool_source_unmatched.tsv",
        # Main pipeline outputs
        fasta    = f"{RESULTS}/output.fasta",
        bed      = f"{RESULTS}/output.bed",
        report   = f"{RESULTS}/report.html",
        unmatched = f"{RESULTS}/pattern_unmatched.tsv",
        not_found = f"{RESULTS}/matched_not_found.tsv",
        unresolv = f"{RESULTS}/unresolved.tsv",
        ambig    = f"{RESULTS}/ambiguous.tsv",
