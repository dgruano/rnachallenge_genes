"""
scripts/resolve_ncbi_assembly_accessions.py
Stage 2c — NCBI Assembly Accession Resolution
==============================================

Enriches NCBI transcripts with NC_/NW_ sequence accessions by:
  1. Mapping NC_/NW_ → parent GCF_/GCA_ assembly accessions
  2. Downloading GTF files for each assembly
  3. Extracting genomic coordinates from GTF
  4. Filling in empty chrom fields and updating assembly_accession

Input
-----
results/ncbi_ucsc_resolved.tsv — NCBI transcripts from resolve_ids
  (may have NC_/NW_ accessions and missing chrom values)

Output
------
results/ncbi_assembly_resolved.tsv — enriched with:
  - assembly_accession: updated from NC_/NW_ to GCF_/GCA_
  - chrom: filled from GTF if was empty
  - start, end: updated from GTF for accuracy
  - is_ambiguous: True if coordinates were ambiguous

results/ncbi_assembly_unresolved.tsv — rows that could not be resolved
"""

import gzip
import re
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Optional

import pandas as pd
from Bio import Entrez, SeqIO

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger
from ncbi_assembly_utils import (
    set_entrez_credentials,
    resolve_assembly_ftp,
    download_gtf,
    extract_all_from_gtf,
    extract_annotations_by_geneid,
)

# ── Snakemake interface ──────────────────────────────────────────────────────
log = get_logger("resolve_ncbi_assembly_accessions", snakemake.log[0])
input_resolved = snakemake.input.resolved
out_resolved = snakemake.output.resolved
out_unresolved = snakemake.output.unresolved
out_ambiguous = snakemake.output.ambiguous
cfg = snakemake.config

set_entrez_credentials(cfg["ncbi_email"], cfg.get("ncbi_api_key"))

CACHE_DIR = Path(cfg["cache_dir"])
MAX_RETRIES = int(cfg.get("max_retries", 3))
RETRY_WAIT = float(cfg.get("retry_wait_seconds", 0.5))
BATCH_SIZE = int(cfg.get("ncbi_batch_size", 50))
EFETCH_BATCH_SIZE = 50

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

# Step 3: Download GTFs for all unique assemblies
unique_asms = df_mapped[df_mapped["assembly_accession"].notna()]["assembly_accession"].unique()
log.info(f"Downloading GTFs for {len(unique_asms)} unique assembly(ies)")

asm_ftp_info = resolve_assembly_ftp(list(unique_asms), log=log, max_retries=MAX_RETRIES, retry_wait=RETRY_WAIT)
log.info(f"Resolved FTP paths for {len(asm_ftp_info)} assembly(ies)")

for asm_acc, ftp_info in asm_ftp_info.items():
    urls = ftp_info["urls"]
    gtf_path = download_gtf(asm_acc, urls, CACHE_DIR, log=log, max_retries=MAX_RETRIES, retry_wait=RETRY_WAIT)
    if gtf_path is None:
        log.error(f"Failed to download GTF for {asm_acc}")

# Step 4: Extract coordinates from GTFs
log.info("Extracting genomic coordinates from GTFs")

# Group rows by assembly for efficient GTF parsing
rows_by_asm = defaultdict(list)
for idx, row in df_mapped.iterrows():
    asm_acc = row["assembly_accession"]
    if pd.notna(asm_acc):
        rows_by_asm[asm_acc].append((idx, row))

# Process each assembly
extracted_coords = {}  # {(idx, transcript_id): {chrom, start, end, strand, gene_id, gene_symbol}}
for asm_acc, rows_list in rows_by_asm.items():
    gtf_path = CACHE_DIR / asm_acc / "genomic.gtf.gz"
    if not gtf_path.exists():
        log.warning(f"GTF not found for {asm_acc}: {gtf_path}")
        continue

    # Collect transcript IDs to search for
    tx_ids = set()
    idx_map = {}  # {transcript_id: [(idx, row), ...]}
    for idx, row in rows_list:
        tx_id = row["transcript_id"]
        tx_ids.add(tx_id)
        if tx_id not in idx_map:
            idx_map[tx_id] = []
        idx_map[tx_id].append((idx, row))

    log.info(f"  Extracting coordinates for {len(tx_ids)} transcript(s) from {asm_acc}…")

    # Try direct transcript ID lookup first
    extracted = extract_all_from_gtf(gtf_path, tx_ids, log=log)
    log.info(f"    Matched {len(extracted)}/{len(tx_ids)} by transcript ID")

    # Record matches
    for tx_id, coords in extracted.items():
        for idx, row in idx_map[tx_id]:
            extracted_coords[(idx, tx_id)] = coords

    # Fallback: for unmatched transcripts, try gene ID matching
    unmatched_tx = tx_ids - set(extracted.keys())
    if unmatched_tx:
        log.debug(f"    {len(unmatched_tx)} transcript(s) not found by ID; trying gene ID fallback…")
        geneid_to_tx = defaultdict(list)
        for tx_id in unmatched_tx:
            # Extract gene_id from input rows
            for idx, row in idx_map[tx_id]:
                if pd.notna(row.get("gene_id")):
                    geneid_to_tx[str(row["gene_id"])].append(tx_id)

        if geneid_to_tx:
            extracted_by_gene = extract_annotations_by_geneid(gtf_path, geneid_to_tx, log=log)
            log.info(f"    Matched {len(extracted_by_gene)} additional transcript(s) by gene ID (ambiguous)")
            for tx_id, coords in extracted_by_gene.items():
                for idx, row in idx_map[tx_id]:
                    extracted_coords[(idx, tx_id)] = {**coords, "is_ambiguous": True}

# Step 5: Update dataframe with extracted coordinates
for idx, row in df_mapped.iterrows():
    tx_id = row["transcript_id"]
    asm_acc = row["assembly_accession"]
    key = (idx, tx_id)

    if key in extracted_coords:
        coords = extracted_coords[key]
        if pd.isna(row["chrom"]) or row["chrom"] == "":
            row["chrom"] = coords.get("chrom", "")
        if pd.isna(row["start"]):
            row["start"] = coords.get("start")
        if pd.isna(row["end"]):
            row["end"] = coords.get("end")
        if row["is_ambiguous"] is None or pd.isna(row["is_ambiguous"]):
            row["is_ambiguous"] = coords.get("is_ambiguous", False)
        df_mapped.loc[idx] = row

# Step 6: Write outputs
log.info(f"Writing {len(df_mapped)} resolved row(s) to {out_resolved}")
df_mapped[RESOLVED_COLS].to_csv(out_resolved, sep="\t", index=False)

df_ambiguous = df_mapped[df_mapped["is_ambiguous"] == True][RESOLVED_COLS]
log.info(f"Writing {len(df_ambiguous)} ambiguous row(s) to {out_ambiguous}")
df_ambiguous.to_csv(out_ambiguous, sep="\t", index=False)

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
