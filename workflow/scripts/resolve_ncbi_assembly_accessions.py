"""
scripts/resolve_ncbi_assembly_accessions.py
Stage 2 NCBI — Assembly Accession and URL Resolution
=====================================================

Enriches NCBI transcripts with NC_/NW_ sequence accessions by:
  1. Mapping NC_/NW_ → parent GCF_/GCA_ assembly accessions
  2. Resolving FTP URLs (fasta_url, gtf_url) for each GCF_/GCA_ assembly
  3. Propagating assembly_name, assembly_accession, fasta_url, gtf_url,
     gtf_format into the per-transcript resolved table

Genomic coordinate extraction (previously Steps 3-5) has been moved to
Stage 5 (extract_sequences.py), which downloads the GTF on demand and
fills missing coords before extraction.

Input
-----
results/ncbi_ucsc_resolved.tsv — NCBI transcripts from resolve_ids
  (may have NC_/NW_ accessions and missing chrom values)

Output
------
results/ncbi_assembly_resolved.tsv — enriched with:
  - assembly_accession: updated from NC_/NW_ to GCF_/GCA_
  - assembly_name: human-readable name (e.g., GRCh38.p14)
  - fasta_url: NCBI FTP FASTA URL
  - gtf_url: NCBI FTP GTF URL
  - gtf_format: "gtf"

results/ncbi_assembly_unresolved.tsv — rows that could not be resolved
"""

import re
import sys
import time
from pathlib import Path
from typing import Optional

import pandas as pd
from Bio import Entrez, SeqIO

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger
from ncbi_assembly_utils import (
    set_entrez_credentials,
    resolve_assembly_ftp,
)

# ── Snakemake interface ──────────────────────────────────────────────────────
log = get_logger("resolve_ncbi_assembly_accessions", snakemake.log[0])
input_resolved = snakemake.input.resolved
out_resolved = snakemake.output.resolved
out_unresolved = snakemake.output.unresolved
cfg = snakemake.config

set_entrez_credentials(cfg["ncbi_email"], cfg.get("ncbi_api_key"))

MAX_RETRIES = int(cfg.get("max_retries", 3))
RETRY_WAIT = float(cfg.get("retry_wait_seconds", 0.5))
EFETCH_BATCH_SIZE = 50

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

# ── Patterns for genomic accessions ──────────────────────────────────────────
_GENOMIC_PREFIXES = ("NC_", "NT_", "NW_", "AC_")
_GENOMIC_RE = re.compile(r'\b((?:NC|NT|NW|AC)_\d+(\.\d+)?)\b')


def _extract_assembly_from_genomic_record(gb_record) -> Optional[str]:
    """Extract assembly accession from a genomic GenBank record."""
    for xref in getattr(gb_record, "dbxrefs", []):
        if xref.startswith("Assembly:"):
            return xref.split(":", 1)[1]
    return None


def map_genomic_to_assembly(accessions: list[str]) -> dict[str, Optional[str]]:
    """
    For each genomic accession (NC_/NW_), fetch its parent assembly (GCF_/GCA_).

    Batches up to EFETCH_BATCH_SIZE accessions per efetch call for efficiency.

    Returns: {genomic_accession: assembly_accession or None}
    """
    results: dict[str, Optional[str]] = {acc: None for acc in accessions}
    total = len(accessions)

    log.info(f"Mapping {total} genomic accession(s) to parent assemblies")

    for batch_start in range(0, len(accessions), EFETCH_BATCH_SIZE):
        batch = accessions[batch_start : batch_start + EFETCH_BATCH_SIZE]
        batch_end = min(batch_start + EFETCH_BATCH_SIZE, len(accessions))

        log.info(f"  Batch {batch_start // EFETCH_BATCH_SIZE + 1}: fetching {len(batch)} accessions…")

        for attempt in range(1, MAX_RETRIES + 1):
            try:
                handle = Entrez.efetch(db="nucleotide", id=",".join(batch), rettype="gb", retmode="text")
                records = SeqIO.parse(handle, "genbank")

                for record in records:
                    acc = record.id.split(".")[0] + ("." + record.id.split(".")[-1] if "." in record.id else "")
                    asm_acc = _extract_assembly_from_genomic_record(record)
                    if asm_acc:
                        results[acc] = asm_acc
                        log.debug(f"    {acc} → {asm_acc}")
                    else:
                        log.debug(f"    {acc} → (no assembly found in record)")

                handle.close()
                break
            except Exception as exc:
                log.warning(f"  Batch attempt {attempt}/{MAX_RETRIES} failed: {exc}")
                if attempt < MAX_RETRIES:
                    time.sleep(RETRY_WAIT * attempt)
                else:
                    log.error(f"  All attempts failed for batch {batch_start // EFETCH_BATCH_SIZE + 1}")

        time.sleep(0.02)  # NCBI rate limiting

    return results


# ── Main processing ────────────────────────────────────────────────────────
log.info("Stage 2c: Resolving NCBI assembly accessions")

df = pd.read_csv(input_resolved, sep="\t", dtype={"chrom": "object"})
log.info(f"Loaded {len(df)} NCBI transcript(s)")

# Ensure new schema columns exist (populated in later tasks; NA until then)
for _col in ("assembly_name", "fasta_url", "gtf_url", "gtf_format"):
    if _col not in df.columns:
        df[_col] = pd.NA

# Identify rows with NC_/NW_ accessions (those that need assembly mapping)
needs_assembly = df[df["assembly_accession"].str.startswith(("NC_", "NW_"), na=False)]
log.info(f"Found {len(needs_assembly)} row(s) with NC_/NW_ sequence accessions needing assembly mapping")

if len(needs_assembly) == 0:
    # No work to do; pass through
    log.info("No NC_/NW_ accessions found; writing through as-is")
    df[RESOLVED_COLS].to_csv(out_resolved, sep="\t", index=False)
    pd.DataFrame(columns=["transcript_id", "db_source", "reason"]).to_csv(
        out_unresolved, sep="\t", index=False
    )
    log.info("Done.")
    sys.exit(0)

# Step 1: Map each unique NC_/NW_ → GCF_/GCA_
unique_genomic = needs_assembly["assembly_accession"].unique()
genomic_to_asm = map_genomic_to_assembly(list(unique_genomic))

mapped_count = sum(1 for v in genomic_to_asm.values() if v is not None)
log.info(f"Successfully mapped {mapped_count}/{len(unique_genomic)} genomic accessions to assemblies")

# Step 2: Update the dataframe with mapped assemblies
resolved_rows = []
unresolved_rows = []

for idx, row in needs_assembly.iterrows():
    genomic_acc = row["assembly_accession"]
    asm_acc = genomic_to_asm.get(genomic_acc)

    if asm_acc is None:
        log.debug(f"  {row['transcript_id']}: assembly mapping failed for {genomic_acc}")
        unresolved_rows.append({
            "transcript_id": row["transcript_id"],
            "db_source": row["db_source"],
            "reason": f"assembly_mapping_failed:{genomic_acc}",
        })
        continue

    row_copy = row.copy()
    row_copy["assembly_accession"] = asm_acc
    resolved_rows.append(row_copy)

df_mapped = pd.concat([df[~df.index.isin(needs_assembly.index)]] + [pd.DataFrame(resolved_rows)], ignore_index=True)
log.info(f"Mapped {len(resolved_rows)} row(s); {len(unresolved_rows)} unresolvable")

if len(df_mapped) == 0:
    log.error("No rows to process after assembly mapping")
    pd.DataFrame(columns=RESOLVED_COLS).to_csv(out_resolved, sep="\t", index=False)
    pd.DataFrame(unresolved_rows, columns=["transcript_id", "db_source", "reason"]).to_csv(
        out_unresolved, sep="\t", index=False
    )
    sys.exit(0)

# Step 3: Resolve FTP URLs (fasta_url, gtf_url) for all unique GCF_/GCA_ assemblies
unique_asms = df_mapped[df_mapped["assembly_accession"].notna()]["assembly_accession"].unique()
log.info(f"Resolving FTP URLs for {len(unique_asms)} unique assembly(ies)")

asm_ftp_info = resolve_assembly_ftp(list(unique_asms), log=log, max_retries=MAX_RETRIES, retry_wait=RETRY_WAIT)
log.info(f"Resolved FTP info for {len(asm_ftp_info)} assembly(ies)")

# Step 4: Propagate assembly_name, fasta_url, gtf_url, gtf_format into df_mapped
for idx, row in df_mapped.iterrows():
    asm_acc = row.get("assembly_accession")
    if pd.isna(asm_acc) or asm_acc not in asm_ftp_info:
        continue
    info = asm_ftp_info[asm_acc]
    df_mapped.at[idx, "gtf_url"] = info["gtf_url"]
    df_mapped.at[idx, "fasta_url"] = info["fasta_url"]
    df_mapped.at[idx, "gtf_format"] = "gtf"
    if pd.isna(row.get("assembly_name")):
        df_mapped.at[idx, "assembly_name"] = info.get("assembly_name", asm_acc)

# Step 5: Write outputs
log.info(f"Writing {len(df_mapped)} resolved row(s) to {out_resolved}")
df_mapped[RESOLVED_COLS].to_csv(out_resolved, sep="\t", index=False)

unresolved_df = pd.DataFrame(unresolved_rows)
if len(unresolved_df) > 0:
    log.warning(f"Writing {len(unresolved_df)} unresolved row(s) to {out_unresolved}")
    unresolved_df.to_csv(out_unresolved, sep="\t", index=False)
else:
    log.info("No unresolved rows")
    pd.DataFrame(columns=["transcript_id", "db_source", "reason"]).to_csv(
        out_unresolved, sep="\t", index=False
    )

log.info("Done.")
