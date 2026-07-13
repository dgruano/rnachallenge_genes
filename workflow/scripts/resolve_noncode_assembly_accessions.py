"""
scripts/resolve_noncode_assembly_accessions.py
Stage 3 - NONCODE Assembly Accession Resolution
================================================

Maps assembly build names (UCSC genome names like tair10, ce10, etc.)
in NONCODE-resolved rows to their corresponding NCBI assembly accessions
(GCF_/GCA_).

Input
-----
results/noncode_resolved.tsv — NONCODE-resolved transcripts with UCSC
                                build names in assembly_accession

Output
------
results/noncode_assembly_resolved.tsv — rows with assembly accessions mapped
results/noncode_assembly_unresolved.tsv — rows where UCSC mapping failed
"""

import sys
from pathlib import Path
from typing import Optional

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake interface ──────────────────────────────────────────────────────
# Only initialize snakemake variables when running under Snakemake
if "snakemake" in globals():
    log = get_logger("resolve_noncode_assembly_accessions", snakemake.log[0])
    input_resolved = snakemake.input.resolved
    out_resolved = snakemake.output.resolved
    out_unresolved = snakemake.output.unresolved

# ── Hardcoded UCSC → GCF_/GCA_ mapping ────────────────────────────────────
# Mapping of UCSC genome build names to NCBI assembly accessions
# Based on the 10 species present in NONCODE dataset
# Verified against NCBI Assembly database: https://www.ncbi.nlm.nih.gov/assembly/
UCSC_TO_GCF_MAPPING = {
    "TAIR10": "GCF_000001735.4",  # Arabidopsis thaliana (TAIR10.1)
    "CE10": "GCF_000002985.6",  # Caenorhabditis elegans (WBcel235)
    "DM6": "GCF_000001215.4",  # Drosophila melanogaster (Release 6 plus ISO1 MT)
    "RN6": "GCF_000001895.5",  # Rattus norvegicus (Rnor_6.0)
    "MONDOM5": "GCF_000002295.2",  # Monodelphis domesticus (MonDom5)
    "PONABE2": "GCF_000001545.5",  # Pongo abelii (P_pygmaeus_2.0.2)
    "GALGAL4": "GCF_000002315.6",  # Gallus gallus (GRCg6a)
    "ORNANA1": "GCF_000002275.2",  # Ornithorhynchus anatinus (ASM227v2)
    "BOSTAU6": "GCF_000003055.6",  # Bos taurus (Bos_taurus_UMD_3.1.1)
    "DANRER10": "GCF_000002035.6",  # Danio rerio (GRCz11)
}

RESOLVED_COLS = [
    "transcript_id",
    "db_source",
    "gene_id",
    "gene_symbol",
    "organism",
    "assembly_name",
    "assembly_accession",
    "fasta_url",
    "gtf_url",
    "gtf_format",
    "chrom",
    "start",
    "end",
    "strand",
    "is_ambiguous",
]


def _validate_ucsc_mapping():
    """Verify no duplicate GCF accessions across species.

    This validation ensures that each UCSC genome build maps to a unique
    NCBI assembly accession, preventing data corruption from incorrect mappings.

    Raises:
        ValueError: If duplicate GCF accessions are found in the mapping.
    """
    accessions = list(UCSC_TO_GCF_MAPPING.values())
    if len(accessions) != len(set(accessions)):
        duplicates = [acc for acc in accessions if accessions.count(acc) > 1]
        raise ValueError(
            f"Duplicate GCF accessions found in UCSC_TO_GCF_MAPPING: {set(duplicates)}"
        )


# Validate mapping at module load time
_validate_ucsc_mapping()


def normalize_ucsc_name(build: str) -> str:
    """Normalize UCSC assembly name to uppercase, remove whitespace."""
    if build is None or pd.isna(build):
        return ""
    return str(build).strip().upper()


def is_ucsc_assembly_accession(accession: str) -> bool:
    """
    Check if accession is a known UCSC genome build name.

    This function checks dictionary membership against UCSC_TO_GCF_MAPPING,
    NOT pattern matching. Only UCSC names explicitly in the mapping are
    recognized as valid UCSC accessions.

    Args:
        accession: String to check, e.g., "ce10", "dm6", "tair10"

    Returns:
        True if accession (normalized) exists in UCSC_TO_GCF_MAPPING, False otherwise.
    """
    if accession is None or pd.isna(accession):
        return False
    normalized = normalize_ucsc_name(accession)
    return normalized in UCSC_TO_GCF_MAPPING


def map_ucsc_to_gcf(assembly_name: str) -> Optional[str]:
    """
    Map UCSC genome build name to GCF_/GCA_ accession.

    Args:
        assembly_name: e.g., "tair10", "ce10", "dm6"

    Returns:
        GCF_/GCA_ accession if mapping exists, None otherwise
    """
    if not assembly_name:
        return None
    normalized = normalize_ucsc_name(assembly_name)
    return UCSC_TO_GCF_MAPPING.get(normalized)


def filter_and_resolve_noncode(
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Filter NONCODE rows with UCSC assembly names and attempt mapping.

    Args:
        df: Input dataframe with noncode_resolved data

    Returns:
        Tuple of (resolved_df, unresolved_df)
    """
    if df.empty:
        return df.copy(), pd.DataFrame()

    resolved_rows = []
    unresolved_rows = []

    for idx, row in df.iterrows():
        # Note: pandas Series.get() returns None if key is missing
        # (unlike dict.get which requires default parameter)
        assembly_val = row.get("assembly_accession")

        # If not a known UCSC name, pass through unchanged (for resolved)
        if not is_ucsc_assembly_accession(assembly_val):
            resolved_rows.append(row.to_dict())
            continue

        # Try to map UCSC to GCF_
        mapped = map_ucsc_to_gcf(str(assembly_val))
        if mapped:
            row_dict = row.to_dict()
            row_dict["assembly_name"] = str(
                assembly_val
            )  # UCSC name as human-readable ID
            row_dict["assembly_accession"] = mapped
            resolved_rows.append(row_dict)
        else:
            # Could not map - goes to unresolved
            unresolved_rows.append(
                {
                    "transcript_id": row.get("transcript_id"),
                    "db_source": row.get("db_source", "noncode"),
                    "reason": f"ucsc_mapping_failed:{assembly_val}",
                }
            )

    resolved = pd.DataFrame(resolved_rows) if resolved_rows else pd.DataFrame()
    unresolved = pd.DataFrame(unresolved_rows) if unresolved_rows else pd.DataFrame()

    return resolved, unresolved


# ── Main processing ────────────────────────────────────────────────────────
# Only run when executed by Snakemake
if "snakemake" in globals():
    log.info("Stage 3: Resolving NONCODE assembly accessions")

    try:
        df = pd.read_csv(input_resolved, sep="\t", dtype={"chrom": "object"})
    except FileNotFoundError as e:
        log.error(f"Input file not found: {input_resolved}")
        raise
    except pd.errors.ParserError as e:
        log.error(f"Failed to parse TSV: {e}")
        raise

    log.info(f"Loaded {len(df)} NONCODE transcript(s)")

    # Ensure new schema columns exist
    for _col in ("assembly_name", "fasta_url", "gtf_url", "gtf_format"):
        if _col not in df.columns:
            df[_col] = pd.NA

    # Validate required columns
    required_columns = ["transcript_id", "assembly_accession"]
    missing = [c for c in required_columns if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Filter and resolve using the dedicated function
    resolved, unresolved = filter_and_resolve_noncode(df)

    # Log statistics
    ucsc_mask = df["assembly_accession"].apply(is_ucsc_assembly_accession)
    needs_mapping = df[ucsc_mask]
    log.info(
        f"Found {len(needs_mapping)} row(s) with UCSC assembly names needing mapping"
    )

    # Write outputs
    log.info(f"Writing {len(resolved)} resolved row(s) to {out_resolved}")
    try:
        if len(resolved) > 0:
            # Ensure columns exist and are in the right order (handle cases with extra columns)
            cols_to_write = [c for c in RESOLVED_COLS if c in resolved.columns]
            # If there are extra columns, keep them at the end
            extra_cols = [c for c in resolved.columns if c not in RESOLVED_COLS]
            resolved[cols_to_write + extra_cols].to_csv(
                out_resolved, sep="\t", index=False
            )
        else:
            pd.DataFrame(columns=RESOLVED_COLS).to_csv(
                out_resolved, sep="\t", index=False
            )
    except IOError as e:
        log.error(f"Failed to write resolved output to {out_resolved}: {e}")
        raise

    try:
        if len(unresolved) > 0:
            log.warning(
                f"Writing {len(unresolved)} unresolved row(s) to {out_unresolved}"
            )
            unresolved.to_csv(out_unresolved, sep="\t", index=False)
        else:
            log.info("No unresolved rows")
            pd.DataFrame(columns=["transcript_id", "db_source", "reason"]).to_csv(
                out_unresolved, sep="\t", index=False
            )
    except IOError as e:
        log.error(f"Failed to write unresolved output to {out_unresolved}: {e}")
        raise

    log.info(f"Resolved: {len(resolved)}, Unresolved: {len(unresolved)}")
    log.info("Done.")
