# ============================================================
# Rules: source_mapping
# ============================================================
# Stage 0 — Download tool FASTAs and map transcript IDs.
#
# Rule graph:
#   download_tool_datasets
#         └── map_ids_to_tools
#                   └── (feeds into parse_ids as side-channel annotation)
# ============================================================

# ── Rule 0a: Download ────────────────────────────────────────
rule download_tool_datasets:
    """
    Download FASTA files from each of the 15 lncRNA/mRNA classifiers listed
    in config/tool_sources.yaml.  Failed or inaccessible URLs are skipped
    gracefully; a JSON manifest records the outcome of every download attempt.
    """
    input:
        tool_sources = "config/tool_sources.yaml",
    output:
        manifest = f"{RESULTS}/tool_datasets_manifest.json",
    params:
        datasets_dir = config.get("tool_datasets_dir", "resources/tool_datasets"),
    log:
        f"{LOGS}/download_tool_datasets.log",
    benchmark:
        f"{BENCHMARKS}/download_tool_datasets.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 120,   # minutes — network-bound
        mem_mb          = 2048,
        cpus_per_task   = 1,
    script:
        "../scripts/download_tool_datasets.py"


# ── Rule 0b: Map ─────────────────────────────────────────────
rule map_ids_to_tools:
    """
    Cross-reference every transcript ID in RNAChallenge.fa against the IDs
    found in each tool's downloaded datasets.

    Outputs:
      tool_source_map.tsv  — transcript_id → tools (comma-sep), primary_tool
      tool_source_stats.tsv — per-tool match counts
    """
    input:
        challenge_fasta = config["input_fastas"][0],
        manifest        = f"{RESULTS}/tool_datasets_manifest.json",
    output:
        tool_map   = f"{RESULTS}/tool_source_map.tsv",
        tool_stats = f"{RESULTS}/tool_source_stats.tsv",
        unmatched  = f"{RESULTS}/tool_source_unmatched.tsv",
    log:
        f"{LOGS}/map_ids_to_tools.log",
    benchmark:
        f"{BENCHMARKS}/map_ids_to_tools.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 30,
        mem_mb          = 4096,
        cpus_per_task   = 1,
    script:
        "../scripts/map_ids_to_tools.py"
