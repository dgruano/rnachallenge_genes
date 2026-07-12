# ============================================================
# Rules: per-accession assembly downloads (fan-out/aggregate)
# ============================================================
# Stage 3 / Phase 4 — Download & Cache NCBI Assemblies
#
# DAG:
#   prepare_accession_list  →  download_assembly (×N)  →  download_assemblies_done
#
# Fault tolerance: download_assembly always exits 0 and writes a status
# sentinel (.download_done = "ok" or "failed: reason"). This means a single
# failed download never blocks the aggregate rule or downstream stages.
# Failed assemblies are reported in unresolved_assemblies.tsv with the
# failure mode extracted from the per-accession log file.
#
# Cache layout:
#   resources/cache/
#     <accession>/
#       .download_done     (Snakemake output — always written)
#       genome.fasta       (created by script on success)
#       genome.fasta.fai   (created by script on success)
#
# Output:
#   results/assembly_accessions.txt   — deduplicated GCF_/GCA_ list
#   results/downloaded_assemblies.tsv — successfully downloaded rows
#   results/unresolved_assemblies.tsv — non-GCF_/GCA_ + failed rows with fail_detail
#   results/.assemblies_ready         — sentinel for extract_sequences
# ============================================================


def get_assemblies(wildcards=None):
    """Return list of GCF_/GCA_ accessions once prepare_accession_list has run."""
    ck = checkpoints.prepare_accession_list.get()
    with open(ck.output.accession_list) as fh:
        return [line.strip() for line in fh if line.strip()]


checkpoint prepare_accession_list:
    input:
        resolved = f"{RESULTS}/ncbi_chromosome_resolved.tsv",
    output:
        accession_list = f"{RESULTS}/assembly_accessions.txt",
    log:
        f"{LOGS}/prepare_accession_list.log",
    benchmark:
        f"{BENCHMARKS}/prepare_accession_list.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 5,
        mem_mb          = 512,
        cpus_per_task   = 1,
    run:
        from pathlib import Path
        import pandas as pd
        df = pd.read_csv(input.resolved, sep="\t")
        accessions = (
            df["assembly_accession"]
            .dropna()
            .astype(str)
            .str.strip()
            .loc[lambda s: s.str.match(r"GC[FA]_\d+\.\d+")]
            .drop_duplicates()
            .tolist()
        )
        Path(output.accession_list).write_text("\n".join(accessions) + "\n")


rule download_assembly:
    wildcard_constraints:
        accession = r"GC[FA]_\d+\.\d+",
    output:
        status = f"{CACHE}/{{accession}}/.download_done",
    log:
        f"{LOGS}/download_assembly/{{accession}}.log",
    benchmark:
        f"{BENCHMARKS}/download_assembly/{{accession}}.tsv",
    resources:
        slurm_partition  = "compute",
        runtime          = 120,
        mem_mb           = 2048,
        cpus_per_task    = 1,
        ncbi_connections = 1,
    envmodules:
        "samtools/1.23.1"
    script:
        "../scripts/download_assembly.py"


rule download_assemblies_done:
    input:
        status_files   = lambda wc: expand(
            f"{CACHE}/{{accession}}/.download_done",
            accession=get_assemblies(wc),
        ),
        log_files      = lambda wc: expand(
            f"{LOGS}/download_assembly/{{accession}}.log",
            accession=get_assemblies(wc),
        ),
        accession_list = f"{RESULTS}/assembly_accessions.txt",
        resolved       = f"{RESULTS}/ncbi_chromosome_resolved.tsv",
    output:
        done       = f"{RESULTS}/.assemblies_ready",
        downloaded = f"{RESULTS}/downloaded_assemblies.tsv",
        unresolved = f"{RESULTS}/unresolved_assemblies.tsv",
    log:
        f"{LOGS}/download_assemblies_done.log",
    benchmark:
        f"{BENCHMARKS}/download_assemblies_done.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 5,
        mem_mb          = 512,
        cpus_per_task   = 1,
    script:
        "../scripts/aggregate_downloads.py"
