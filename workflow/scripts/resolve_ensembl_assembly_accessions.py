"""
scripts/resolve_ensembl_assembly_accessions.py
Stage 2d - Ensembl Assembly Accession Resolution
==================================================

Maps assembly build names (GRCh38, GRCz11, etc.) in Ensembl-resolved rows
to their corresponding NCBI assembly accessions (GCF_/GCA_).

Input
-----
results/ensembl_resolved.tsv — Ensembl-resolved transcripts with possible
                                GRC* build names in assembly_accession

Output
------
results/ensembl_assembly_resolved.tsv — rows with assembly accessions mapped
                                        or pass-through if not GRC*
results/ensembl_assembly_unresolved.tsv — rows where GRC* mapping failed
"""

import sys
from pathlib import Path
from typing import Optional

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake interface ──────────────────────────────────────────────────────
log = get_logger("resolve_ensembl_assembly_accessions", snakemake.log[0])
input_resolved = snakemake.input.resolved
out_resolved = snakemake.output.resolved
out_unresolved = snakemake.output.unresolved

# ── Hardcoded GRC* → GCF_/GCA_ mapping ───────────────────────────────────────
ASSEMBLY_NAME_MAPPING = {
    "GRCH38": "GCF_000001405.40",    # Homo sapiens (human)
    "GRCH37": "GCF_000001405.39",    # Homo sapiens GRCh37 (older)
    "GRCZ11": "GCF_000002035.6",     # Danio rerio (zebrafish)
    "GRCZ10": "GCF_000002035.5",     # Danio rerio (older)
    "GRCRH1": "GCF_000008735.2",     # Macaca mulatta (rhesus macaque)
    "GRCM39": "GCF_000001635.27",    # Mus musculus (mouse)
    "GRCM38": "GCF_000001635.26",    # Mus musculus (older)
    "GRCH13": "GCF_000004545.3",     # Ciona intestinalis (sea squirt)
    "BDGP6": "GCF_000001215.4",      # Drosophila melanogaster (fruit fly)
    "BDGP5": "GCF_000001215.3",      # Drosophila melanogaster (older)
    "BDGP6_32": "GCF_000001215.4",   # Drosophila release variant
}

RESOLVED_COLS = [
    "transcript_id",
    "db_source",
    "gene_id",
    "gene_symbol",
    "organism",
    "assembly_accession",
    "chrom",
    "start",
    "end",
    "strand",
    "is_ambiguous",
]


def normalize_build_name(build: str) -> str:
    """Normalize assembly name to uppercase, remove whitespace."""
    if build is None or pd.isna(build):
        return ""
    return str(build).strip().upper()


def is_grc_assembly_accession(accession: str) -> bool:
    """Check if accession matches GRC* pattern."""
    if accession is None or pd.isna(accession):
        return False
    normalized = str(accession).strip().upper()
    return normalized.startswith("GRC")


def map_grc_to_gcf(assembly_name: str) -> Optional[str]:
    """
    Map GRC* assembly name to GCF_/GCA_ accession.

    Args:
        assembly_name: e.g., "GRCh38", "GRCz11"

    Returns:
        GCF_/GCA_ accession if mapping exists, None otherwise
    """
    if not assembly_name:
        return None
    normalized = normalize_build_name(assembly_name)
    return ASSEMBLY_NAME_MAPPING.get(normalized)


# ── Main processing ────────────────────────────────────────────────────────
log.info("Stage 2d: Resolving Ensembl assembly accessions")

df = pd.read_csv(input_resolved, sep="\t", dtype={"chrom": "object"})
log.info(f"Loaded {len(df)} Ensembl transcript(s)")

# Identify rows with GRC* accessions (those that need mapping)
grc_mask = df["assembly_accession"].apply(is_grc_assembly_accession)
needs_mapping = df[grc_mask]
log.info(f"Found {len(needs_mapping)} row(s) with GRC* assembly names needing mapping")

# Process all rows
resolved_rows = []
unresolved_rows = []

for idx, row in df.iterrows():
    assembly_val = row.get("assembly_accession")

    # If not a GRC* accession, pass through unchanged (for resolved)
    if not is_grc_assembly_accession(assembly_val):
        resolved_rows.append(row.to_dict())
        continue

    # Try to map GRC* to GCF_
    mapped = map_grc_to_gcf(str(assembly_val))
    if mapped:
        row_dict = row.to_dict()
        row_dict["assembly_accession"] = mapped
        log.debug(f"  {row.get('transcript_id')}: {assembly_val} → {mapped}")
        resolved_rows.append(row_dict)
    else:
        # Could not map - goes to unresolved
        log.warning(f"  {row.get('transcript_id')}: unmappable GRC* accession: {assembly_val}")
        unresolved_rows.append({
            "transcript_id": row.get("transcript_id"),
            "db_source": row.get("db_source", "ensembl"),
            "reason": f"grc_mapping_failed:{assembly_val}",
        })

# Construct output dataframes
resolved = pd.DataFrame(resolved_rows) if resolved_rows else pd.DataFrame()
unresolved = pd.DataFrame(unresolved_rows) if unresolved_rows else pd.DataFrame()

# Write outputs
log.info(f"Writing {len(resolved)} resolved row(s) to {out_resolved}")
if len(resolved) > 0:
    # Ensure columns exist and are in the right order (handle cases with extra columns)
    cols_to_write = [c for c in RESOLVED_COLS if c in resolved.columns]
    # If there are extra columns, keep them at the end
    extra_cols = [c for c in resolved.columns if c not in RESOLVED_COLS]
    resolved[cols_to_write + extra_cols].to_csv(out_resolved, sep="\t", index=False)
else:
    pd.DataFrame(columns=RESOLVED_COLS).to_csv(out_resolved, sep="\t", index=False)

if len(unresolved) > 0:
    log.warning(f"Writing {len(unresolved)} unresolved row(s) to {out_unresolved}")
    unresolved.to_csv(out_unresolved, sep="\t", index=False)
else:
    log.info("No unresolved rows")
    pd.DataFrame(columns=["transcript_id", "db_source", "reason"]).to_csv(
        out_unresolved, sep="\t", index=False
    )

log.info(f"Resolved: {len(resolved)}, Unresolved: {len(unresolved)}")
log.info("Done.")