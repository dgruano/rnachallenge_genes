"""
scripts/resolve_ncbi_chromosome_accessions.py
Post-merge resolver: NC_/NT_/NW_ → parent GCF_ assembly
=========================================================

Reads resolved_ids.tsv (output of merge_resolved) and maps any remaining
chromosomal accessions (NC_/NT_/NW_) to their parent GCF_/GCA_ assembly
accession via NCBI elink.

This is a catch-all pass between merge_resolved and download_assemblies:
any db_source that contributed chromosomal accessions gets fixed here.

Input
-----
results/resolved_ids.tsv

Output
------
results/ncbi_chromosome_resolved.tsv — full TSV with NC_/NT_/NW_ rows patched
results/ncbi_chromosome_unresolved.tsv — rows whose chromosomal accession
                                          could not be mapped to a parent assembly
"""

import sys
import time
from pathlib import Path
from typing import Optional

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger
from ncbi_assembly_utils import set_entrez_credentials, map_genomic_to_assembly_elink

_CHROMOSOMAL_PREFIXES = ("NC_", "NT_", "NW_", "AC_")


def is_chromosomal_accession(value) -> bool:
    """Return True if value is a chromosomal-level NCBI accession (NC_/NT_/NW_/AC_)."""
    try:
        if value is None or pd.isna(value):
            return False
    except (TypeError, ValueError):
        pass
    s = str(value).strip()
    return s.startswith(_CHROMOSOMAL_PREFIXES)


def resolve_chromosomal_rows(
    df: pd.DataFrame,
    mapping: dict[str, Optional[str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split df into (resolved, unresolved) using the pre-computed mapping.

    Rows whose assembly_accession is already GCF_/GCA_ pass through into
    resolved unchanged. Rows with chromosomal accessions are replaced with
    the mapped GCF_ value; those that failed mapping go to unresolved.

    Parameters
    ----------
    df : pd.DataFrame
        Full resolved_ids.tsv content
    mapping : dict
        {chromosomal_accession: gca_accession_or_None}

    Returns
    -------
    (resolved_df, unresolved_df)
        unresolved_df has columns: transcript_id, db_source, reason
    """
    if df.empty:
        return df.copy(), pd.DataFrame(columns=["transcript_id", "db_source", "reason"])

    resolved_rows = []
    unresolved_rows = []

    for _, row in df.iterrows():
        acc = row.get("assembly_accession")
        if not is_chromosomal_accession(acc):
            resolved_rows.append(row.to_dict())
            continue

        gcf = mapping.get(str(acc).strip())
        if gcf:
            d = row.to_dict()
            d["assembly_accession"] = gcf
            resolved_rows.append(d)
        else:
            unresolved_rows.append({
                "transcript_id": row.get("transcript_id"),
                "db_source": row.get("db_source", ""),
                "reason": f"chromosomal_mapping_failed:{acc}",
            })

    resolved = pd.DataFrame(resolved_rows) if resolved_rows else pd.DataFrame(columns=df.columns)
    unresolved = (
        pd.DataFrame(unresolved_rows)
        if unresolved_rows
        else pd.DataFrame(columns=["transcript_id", "db_source", "reason"])
    )
    return resolved, unresolved


# ── Snakemake entry point ────────────────────────────────────────────────────
if "snakemake" in globals():
    log = get_logger("resolve_ncbi_chromosome_accessions", snakemake.log[0])
    input_tsv = snakemake.input.resolved
    out_resolved = snakemake.output.resolved
    out_unresolved = snakemake.output.unresolved
    cfg = snakemake.config

    set_entrez_credentials(cfg["ncbi_email"], cfg.get("ncbi_api_key"))
    MAX_RETRIES = int(cfg.get("max_retries", 3))
    RETRY_WAIT = float(cfg.get("retry_wait_seconds", 0.5))

    log.info("Post-merge resolver: NC_/NT_/NW_ → parent GCF_ assembly")

    df = pd.read_csv(input_tsv, sep="\t", dtype={"chrom": "object"})
    log.info(f"Loaded {len(df)} rows from {input_tsv}")

    # Identify unique chromosomal accessions
    chromosomal_mask = df["assembly_accession"].apply(is_chromosomal_accession)
    unique_chromosomal = (
        df.loc[chromosomal_mask, "assembly_accession"]
        .dropna()
        .unique()
        .tolist()
    )
    log.info(f"Found {chromosomal_mask.sum()} rows with chromosomal accessions "
             f"({len(unique_chromosomal)} unique)")

    if unique_chromosomal:
        mapping = map_genomic_to_assembly_elink(
            unique_chromosomal,
            log=log,
            max_retries=MAX_RETRIES,
            retry_wait=RETRY_WAIT,
        )
        mapped_count = sum(1 for v in mapping.values() if v is not None)
        log.info(f"Mapped {mapped_count}/{len(unique_chromosomal)} chromosomal accessions")
    else:
        log.info("No chromosomal accessions — passing through unchanged")
        mapping = {}

    resolved_df, unresolved_df = resolve_chromosomal_rows(df, mapping)

    log.info(f"Writing {len(resolved_df)} resolved rows to {out_resolved}")
    resolved_df.to_csv(out_resolved, sep="\t", index=False)

    log.info(f"Writing {len(unresolved_df)} unresolved rows to {out_unresolved}")
    unresolved_df.to_csv(out_unresolved, sep="\t", index=False)

    log.info("=" * 60)
    log.info(f"Total input rows            : {len(df)}")
    log.info(f"Chromosomal accessions found: {chromosomal_mask.sum()}")
    log.info(f"Successfully remapped       : {chromosomal_mask.sum() - len(unresolved_df)}")
    log.info(f"Could not remap             : {len(unresolved_df)}")
    log.info("resolve_ncbi_chromosome_accessions complete.")
