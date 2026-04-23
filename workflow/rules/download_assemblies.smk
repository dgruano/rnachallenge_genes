# ============================================================
# Rule: download_assemblies  (checkpoint)
# ============================================================
# Stage 3 / Phase 4 — Download & Cache NCBI Assemblies (Simplified)
#
# Reads the resolved transcript TSV to identify unique assemblies.
# For each GCF_/GCA_ accession:
#   1. Checks if already cached (skips if present)
#   2. Downloads from NCBI FTP → resources/cache/<accession>/
#   3. Indexes with samtools
#
# Non-GCF_/GCA_ accessions are marked as unresolved.
#
# Cache layout:
#   resources/cache/
#     <assembly_accession>/
#       genomic.fna.gz          # Downloaded NCBI FASTA
#       genome.fasta            # Decompressed FASTA
#       genome.fasta.fai        # Index (via samtools faidx)
#
# Output:
#   results/downloaded_assemblies.tsv - successfully downloaded
#   results/unresolved_assemblies.tsv - non-GCF_/GCA_ accessions
#
# A checkpoint is used to allow dependent rules to fan out based
# on which assemblies were downloaded.
# ============================================================

checkpoint download_assemblies:
    input:
        resolved = f"{RESULTS}/resolved_ids.tsv",
    output:
        done = f"{RESULTS}/.assemblies_ready",
        downloaded = f"{RESULTS}/downloaded_assemblies.tsv",
        unresolved = f"{RESULTS}/unresolved_assemblies.tsv",
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
