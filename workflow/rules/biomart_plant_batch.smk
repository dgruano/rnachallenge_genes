# ============================================================
# Rule: biomart_plant_batch
# ============================================================
# Batch query Ensembl Plants BioMart for unresolved plant IDs.
# Much more efficient than REST API for bulk queries.
# ============================================================

rule biomart_plant_batch:
    input:
        unresolved = f"{RESULTS}/external_unresolved.tsv",
    output:
        resolved   = f"{RESULTS}/biomart_resolved.tsv",
        unresolved = f"{RESULTS}/biomart_unresolved.tsv",
    log:
        f"{LOGS}/biomart_plant_batch.log",
    benchmark:
        f"{BENCHMARKS}/biomart_plant_batch.tsv",
    resources:
        slurm_partition = "compute",
        runtime         = 5,
        mem_mb          = 512,
        cpus_per_task   = 1,
    # ponytail: BioMart detached — Ensembl Plants BioMart is unreliable and the
    # plant GTF/Phytozome/Gramene fallbacks cover these IDs off external_*.
    # Emit header-only outputs to keep the DAG satisfied; pass IDs through as
    # unresolved. Re-enable by restoring `script: ../scripts/biomart_plant_batch.py`.
    run:
        import pandas as pd

        resolved_cols = [
            "transcript_id", "db_source", "gene_id", "gene_symbol", "organism",
            "assembly_accession", "assembly_name", "chrom", "start", "end",
            "strand", "ensembl_plants_release", "is_ambiguous",
        ]
        unresolved_cols = ["transcript_id", "raw_header", "source_file", "reason"]

        pd.DataFrame(columns=resolved_cols).to_csv(
            output.resolved, sep="\t", index=False
        )

        df = pd.read_csv(input.unresolved, sep="\t")
        out = pd.DataFrame({
            "transcript_id": df.get("transcript_id", pd.Series(dtype=str)),
            "raw_header": df.get("raw_header", ""),
            "source_file": df.get("source_file", ""),
            "reason": "biomart_detached",
        })
        out.to_csv(output.unresolved, sep="\t", index=False)
