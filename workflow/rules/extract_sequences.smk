# ============================================================
# Rule: extract_sequences
# ============================================================
# Stage 4 — Extract Gene + Flanking Sequences
#
# For each resolved transcript:
#   1. Extends gene coordinates by upstream_bp / downstream_bp
#      (strand-aware, clamped to chromosome bounds)
#   2. Extracts the sequence from the cached indexed assembly
#      using samtools faidx
#   3. Reverse-complements if gene is on minus strand
#   4. Writes FASTA and BED outputs
#
# FASTA header format:
#   >{transcript_id} gene={gene_id} symbol={gene_symbol} organism={organism}
#   assembly={assembly} loc={chrom}:{ext_start}-{ext_end}({strand})
#   upstream={upstream_bp} downstream={downstream_bp}
# ============================================================

rule extract_sequences:
    input:
        resolved = f"{RESULTS}/resolved_ids.tsv",
        sentinel = f"{RESULTS}/.assemblies_ready",
    output:
        fasta = f"{RESULTS}/output.fasta",
        bed   = f"{RESULTS}/output.bed",
        failed  = f"{RESULTS}/extraction_failed.tsv",
    log:
        f"{LOGS}/extract_sequences.log",
    benchmark:
        f"{BENCHMARKS}/extract_sequences.tsv",
    params:
        upstream   = UPSTREAM,
        downstream = DOWNSTREAM,
        cache_dir  = CACHE,
    resources:
        slurm_partition = "compute",
        runtime         = 120,
        mem_mb          = 8192,
        cpus_per_task   = 4,
    script:
        "../scripts/extract_sequences.py"
