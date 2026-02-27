# ============================================================
# Rule: download_assemblies  (checkpoint)
# ============================================================
# Stage 3 — Download & Cache Genome Assemblies
#
# Reads the resolved transcript TSV to identify all unique
# (organism, assembly_accession) pairs. For each pair:
#   1. Checks if already cached (skips if present)
#   2. Determines FTP source (NCBI FTP for GCF_/GCA_, Ensembl FTP otherwise)
#   3. Downloads primary genome FASTA (chromosomes / top-level)
#   4. Decompresses (.gz)
#   5. Indexes with samtools faidx
#
# Cache layout:
#   resources/cache/
#     <assembly_accession>/
#       genome.fasta
#       genome.fasta.fai
#
# A checkpoint is used to allow dependent rules to fan out based
# on which assemblies were downloaded.
# ============================================================

checkpoint download_assemblies:
    input:
        resolved = f"{RESULTS}/resolved_ids.tsv",
    output:
        done = f"{RESULTS}/.assemblies_ready",
    log:
        f"{LOGS}/download_assemblies.log",
    benchmark:
        f"{BENCHMARKS}/download_assemblies.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 480,  # 8 hours for network downloads
        mem_mb          = 4096,
        cpus_per_task   = 2,
    script:
        "../scripts/download_assemblies.py"
