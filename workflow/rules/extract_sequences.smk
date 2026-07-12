# ============================================================
# Rules: per-batch sequence extraction (fan-out/aggregate)
# ============================================================
# Stage 4 — Extract Gene + Flanking Sequences
#
# DAG:
#   split_batches  →  extract_sequences_batch (×N)  →  extract_sequences
#
# Per-batch outputs (results/sequences/):
#   batch_NNNN.fasta
#   batch_NNNN.bed
#   batch_NNNN.failed.tsv
#
# Final outputs:
#   results/output.fasta
#   results/output.bed
#   results/extraction_failed.tsv
# ============================================================


def get_batch_ids(wildcards=None):
    """Return list of batch IDs once split_batches has run."""
    ck = checkpoints.split_batches.get()
    with open(ck.output.manifest) as fh:
        return [line.strip() for line in fh if line.strip()]


checkpoint split_batches:
    input:
        resolved = f"{RESULTS}/ncbi_chromosome_resolved.tsv",
    output:
        manifest  = f"{RESULTS}/batch_manifest.txt",
        batch_dir = directory(f"{RESULTS}/batches"),
    log:
        f"{LOGS}/split_batches.log",
    benchmark:
        f"{BENCHMARKS}/split_batches.tsv",
    params:
        batch_size = config.get("extraction_batch_size", 5000),
    resources:
        slurm_partition = "compute",
        runtime         = 5,
        mem_mb          = 1024,
        cpus_per_task   = 1,
    script:
        "../scripts/split_batches.py"


rule extract_sequences_batch:
    input:
        resolved = f"{RESULTS}/batches/{{batch_id}}.tsv",
        sentinel = f"{RESULTS}/.assemblies_ready",
    output:
        fasta  = f"{RESULTS}/sequences/{{batch_id}}.fasta",
        bed    = f"{RESULTS}/sequences/{{batch_id}}.bed",
        failed = f"{RESULTS}/sequences/{{batch_id}}.failed.tsv",
    log:
        f"{LOGS}/extract_sequences/{{batch_id}}.log",
    benchmark:
        f"{BENCHMARKS}/extract_sequences/{{batch_id}}.tsv",
    params:
        upstream   = UPSTREAM,
        downstream = DOWNSTREAM,
        cache_dir  = CACHE,
    resources:
        slurm_partition = "compute",
        runtime         = 60,
        mem_mb          = 2048,
        cpus_per_task   = 1,
    script:
        "../scripts/extract_sequences.py"


rule extract_sequences:
    input:
        fastas  = lambda wc: expand(
            f"{RESULTS}/sequences/{{batch_id}}.fasta",
            batch_id=get_batch_ids(wc),
        ),
        beds    = lambda wc: expand(
            f"{RESULTS}/sequences/{{batch_id}}.bed",
            batch_id=get_batch_ids(wc),
        ),
        faileds = lambda wc: expand(
            f"{RESULTS}/sequences/{{batch_id}}.failed.tsv",
            batch_id=get_batch_ids(wc),
        ),
    output:
        fasta  = f"{RESULTS}/output.fasta",
        bed    = f"{RESULTS}/output.bed",
        failed = f"{RESULTS}/extraction_failed.tsv",
    log:
        f"{LOGS}/extract_sequences.log",
    benchmark:
        f"{BENCHMARKS}/extract_sequences.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 15,
        mem_mb          = 2048,
        cpus_per_task   = 1,
    script:
        "../scripts/aggregate_sequences.py"
