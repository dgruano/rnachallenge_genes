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
#   results/assembly_download_manifest.tsv — cache_key/fasta_url manifest
#   results/downloaded_assemblies.tsv — successfully downloaded rows
#   results/unresolved_assemblies.tsv — non-GCF_/GCA_ + failed rows with fail_detail
#   results/.assemblies_ready         — sentinel for extract_sequences
# ============================================================


def get_assemblies(wildcards=None):
    """Return list of manifest cache keys once prepare_accession_list has run."""
    ck = checkpoints.prepare_accession_list.get()
    import pandas as pd

    manifest = pd.read_csv(ck.output.accession_list, sep="\t")
    if manifest.empty or "cache_key" not in manifest.columns:
        return []
    return (
        manifest["cache_key"]
        .dropna()
        .astype(str)
        .str.strip()
        .loc[lambda s: s.ne("")]
        .drop_duplicates()
        .tolist()
    )


checkpoint prepare_accession_list:
    input:
        resolved = f"{RESULTS}/ncbi_chromosome_resolved.tsv",
    output:
        accession_list = f"{RESULTS}/assembly_download_manifest.tsv",
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
        import sys
        from pathlib import Path
        import pandas as pd

        sys.path.insert(0, str(Path("workflow/scripts").resolve()))
        from download_manifest_utils import build_download_manifest

        df = pd.read_csv(input.resolved, sep="\t")
        manifest = build_download_manifest(df)
        manifest.to_csv(output.accession_list, sep="\t", index=False)


rule download_assembly:
    input:
        manifest = f"{RESULTS}/assembly_download_manifest.tsv",
    wildcard_constraints:
        accession = r"[A-Za-z0-9._-]+",
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
        # Phytozome genome FASTAs come from a separate JGI-authed fan-out
        # (download_phytozome_fasta.smk); gate the sentinel on them too so
        # extract_sequences waits for them.
        phytozome      = lambda wc: expand(
            f"{CACHE}/phytozome_{{species}}/.download_done",
            species=get_phytozome_species(wc),
        ),
        accession_list = f"{RESULTS}/assembly_download_manifest.tsv",
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
