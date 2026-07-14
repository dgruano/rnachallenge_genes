# ============================================================
# Rule: download_phytozome_fasta (per-species fan-out)
# ============================================================
# Caches the Phytozome genome FASTA (from JGI, same version as the resolved
# GFF3) so extract_sequences can slice phytozome-resolved coordinates.
#
# DAG:
#   resolve_phytozome_gtf → merge → resolved_ids.tsv
#     └─ get_phytozome_species (gated on prepare_accession_list checkpoint)
#          → download_phytozome_fasta (×species)
#              → resources/cache/phytozome_<species>/genome.fasta{,.fai}
#              → (input to) download_assemblies_done → .assemblies_ready
#
# Only species that actually resolved get a FASTA (post-merge, so coordinate-less
# rows dropped by merge are excluded). Fault-tolerant like download_assembly:
# the script always writes a .download_done sentinel and never blocks the run.
#
# NOTE: the fan-out reads results/resolved_ids.tsv (the merged resolved table).
# The download stage otherwise keys off ncbi_chromosome_resolved.tsv, which is a
# misnomer — it's the full merged table with NCBI chrom names patched (rename
# tracked in PIPELINE_AUDIT.md → Refactors). Both carry the same phytozome
# species set; resolved_ids.tsv is the honest name.
# ============================================================


def get_phytozome_species(wildcards=None):
    """Distinct phytozome species (config keys) present in the resolved table."""
    # Gate on the checkpoint so resolved_ids.tsv is final before we read it.
    checkpoints.prepare_accession_list.get()
    import pandas as pd

    df = pd.read_csv(f"{RESULTS}/resolved_ids.tsv", sep="\t")
    if df.empty or "db_source" not in df.columns or "organism" not in df.columns:
        return []
    phyto = df.loc[df["db_source"].astype(str) == "phytozome", "organism"]
    return (
        phyto.dropna()
        .astype(str)
        .str.strip()
        .loc[lambda s: s.ne("")]
        .drop_duplicates()
        .tolist()
    )


# phytozome_<species> keys also match download_assembly's {accession}; this rule
# owns them.
ruleorder: download_phytozome_fasta > download_assembly


rule download_phytozome_fasta:
    input:
        resolved = f"{RESULTS}/resolved_ids.tsv",
    wildcard_constraints:
        species = r"[^/]+",
    output:
        status = f"{CACHE}/phytozome_{{species}}/.download_done",
    log:
        f"{LOGS}/download_phytozome_fasta/{{species}}.log",
    benchmark:
        f"{BENCHMARKS}/download_phytozome_fasta/{{species}}.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 240,
        mem_mb          = 4096,
        cpus_per_task   = 1,
    envmodules:
        "samtools/1.23.1"
    script:
        "../scripts/download_phytozome_fasta.py"
