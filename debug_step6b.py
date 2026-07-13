#!/usr/bin/env python3
"""
debug_step6b.py
===============
Standalone debugger for Step 6b of resolve_abandoned_accessions.py.

Reads results/abandoned_resolved_debug_bak.tsv, extracts unresolved
transcript IDs together with their already-known Gene IDs (from genomic
record tracking), then runs:

  Phase A — batch elink gene→assembly   (batch_link_genes_to_assemblies)
  Phase B — batch esummary UIDs→metadata (resolve_assembly_uids_map)
  Phase C — per-gene nuccore fallback for genes with no elink result,
             using the genomic accession already stored in the TSV
             (fetch_assembly_from_nuccore)

Intentionally skips batch_fetch_gene_info: we already have gene IDs from
Step 2 genomic records, so downloading full Entrezgene XML for 1800+ genes
is wasted work and is the likely cause of the hang in the full pipeline.

Usage
-----
  python debug_step6b.py [--input results/abandoned_resolved_debug_bak.tsv]
                         [--output debug_step6b_results.tsv]
                         [--limit N]   # process only first N rows (testing)
"""

import argparse
import sys
import time
from pathlib import Path

import pandas as pd
from Bio import Entrez

sys.path.insert(0, str(Path(__file__).parent / "workflow" / "scripts"))
from ncbi_entrez_utils import (
    batch_link_genes_to_assemblies,
    fetch_assembly_from_nuccore,
    resolve_assembly_uids_map,
)

# ── Config ────────────────────────────────────────────────────────────────────
NCBI_EMAIL = "daniel.garciaruano@ibgc.cnrs.fr"
NCBI_API_KEY = "fdc3c3e0bf43bfebc561e789a0884b879308"
DELAY = 0.0  # batched requests take >> 0.1s; no artificial sleep needed

DEFAULT_INPUT = "results/abandoned_resolved_debug_bak.tsv"
DEFAULT_OUTPUT = "debug_step6b_results.tsv"

# ── CLI ───────────────────────────────────────────────────────────────────────


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", default=DEFAULT_INPUT)
    p.add_argument("--output", default=DEFAULT_OUTPUT)
    p.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process only the first N unresolved rows (for testing)",
    )
    return p.parse_args()


# ── Helpers ───────────────────────────────────────────────────────────────────


def _parse_geneids(cell) -> list[str]:
    """Split a semicolon-delimited Gene ID cell into a list of strings."""
    if not cell or pd.isna(cell):
        return []
    return [g.strip() for g in str(cell).split(";") if g.strip()]


def _scaffold_from_genomic_acc(genomic_acc: str) -> str:
    """Return the accession if it looks like a scaffold, else empty string."""
    if not genomic_acc or pd.isna(genomic_acc):
        return ""
    _SCAFFOLD_PREFIXES = ("NW_", "NC_", "NT_", "NZ_")
    return (
        genomic_acc.strip() if str(genomic_acc).startswith(_SCAFFOLD_PREFIXES) else ""
    )


# ── Main ──────────────────────────────────────────────────────────────────────


def main():
    args = parse_args()

    Entrez.email = NCBI_EMAIL
    Entrez.api_key = NCBI_API_KEY

    # ── Load debug TSV ────────────────────────────────────────────────────────
    print(f"Reading {args.input} …")
    df = pd.read_csv(args.input, sep="\t", dtype=str)
    print(f"  Total rows: {len(df)}")

    # Filter to rows that need Step 6b: unresolved with no assembly found
    no_asm = df[
        (df["resolution_status"] == "unresolved")
        & (df["unresolved_reason"] == "no_assembly_found")
    ].copy()

    if args.limit:
        no_asm = no_asm.head(args.limit)
        print(f"  (Limited to first {args.limit} rows)")

    print(f"  Rows needing Step 6b: {len(no_asm)}")

    # ── Build gene_id → [transcript_ids] map ─────────────────────────────────
    # Prefer genomic_geneids (from scaffold record, most reliable).
    # Fall back to transcript_geneids if the column exists.
    tx_to_gids: dict[str, list[str]] = {}
    gid_to_txs: dict[str, list[str]] = {}
    tx_to_scaffold: dict[str, str] = {}

    has_tx_geneids = "transcript_geneids" in df.columns

    for _, row in no_asm.iterrows():
        acc = row["transcript_id"]
        gids = _parse_geneids(row.get("genomic_geneids", ""))
        if not gids and has_tx_geneids:
            gids = _parse_geneids(row.get("transcript_geneids", ""))
        tx_to_gids[acc] = gids

        scaffold = _scaffold_from_genomic_acc(row.get("genomic_accession", ""))
        if scaffold:
            tx_to_scaffold[acc] = scaffold

        for gid in gids:
            gid_to_txs.setdefault(gid, []).append(acc)

    all_gene_ids = list(gid_to_txs.keys())
    no_geneid_txs = [acc for acc, gids in tx_to_gids.items() if not gids]

    print(f"\nGene ID summary:")
    print(f"  Transcripts with ≥1 Gene ID : {len(tx_to_gids) - len(no_geneid_txs)}")
    print(f"  Transcripts with NO Gene ID : {len(no_geneid_txs)}")
    print(f"  Unique Gene IDs to query    : {len(all_gene_ids)}")
    print(f"  Transcripts with scaffold   : {len(tx_to_scaffold)}")

    # ── Phase A: elink gene → assembly UIDs ──────────────────────────────────
    print(f"\nPhase A: batch elink gene→assembly ({len(all_gene_ids)} gene IDs) …")
    t0 = time.time()
    gene_asm_link = batch_link_genes_to_assemblies(all_gene_ids, DELAY)
    print(f"  Done in {time.time()-t0:.1f}s")

    linked = {gid: uids for gid, uids in gene_asm_link.items() if uids}
    unlinked = [gid for gid in all_gene_ids if not gene_asm_link.get(gid)]
    print(f"  Genes with assembly link : {len(linked)}/{len(all_gene_ids)}")
    print(f"  Genes with NO link       : {len(unlinked)}")

    # ── Phase B: resolve assembly UIDs → accessions ──────────────────────────
    all_uids = [uid for uids in linked.values() for uid in uids]
    print(f"\nPhase B: batch esummary for {len(set(all_uids))} unique assembly UIDs …")
    t0 = time.time()
    asm_uid_map = resolve_assembly_uids_map(all_uids, DELAY)
    print(f"  Done in {time.time()-t0:.1f}s")
    print(f"  Assembly records resolved: {len(asm_uid_map)}")

    # Build gene_id → assembly_accession from elink results
    gene_to_asm: dict[str, str] = {}
    for gid, uids in linked.items():
        asm_info = asm_uid_map.get(uids[0], {})
        asm_acc = asm_info.get("assembly_accession", "")
        if asm_acc and asm_acc != "N/A":
            gene_to_asm[gid] = asm_acc

    # ── Phase C: nuccore scaffold fallback for unlinked genes ─────────────────
    # Use genomic accession stored in the TSV (already have it — no extra fetch)
    # Find which transcripts still need an assembly and have a scaffold acc
    tx_needing_fallback = []
    for acc in no_asm["transcript_id"]:
        gids = tx_to_gids.get(acc, [])
        if any(gid in gene_to_asm for gid in gids):
            continue  # already resolved via elink
        if acc in tx_to_scaffold:
            tx_needing_fallback.append(acc)

    scaffold_to_asm: dict[str, str] = {}
    if tx_needing_fallback:
        unique_scaffolds = list(
            dict.fromkeys(tx_to_scaffold[a] for a in tx_needing_fallback)
        )
        print(
            f"\nPhase C: nuccore scaffold fallback for {len(unique_scaffolds)} unique scaffold(s) …"
        )
        t0 = time.time()
        for scaffold in unique_scaffolds:
            asm_data = fetch_assembly_from_nuccore(scaffold, DELAY)
            asm_acc = asm_data.get("assembly_accession", "")
            if asm_acc and asm_acc not in ("N/A", ""):
                scaffold_to_asm[scaffold] = asm_acc
                print(f"  {scaffold} → {asm_acc}")
            else:
                print(f"  {scaffold} → NOT FOUND  {asm_data.get('error','')}")
        print(
            f"  Done in {time.time()-t0:.1f}s  ({len(scaffold_to_asm)}/{len(unique_scaffolds)} resolved)"
        )

    # ── Assemble output rows ──────────────────────────────────────────────────
    rows = []
    for _, row in no_asm.iterrows():
        acc = row["transcript_id"]
        gids = tx_to_gids.get(acc, [])

        # Try elink path first
        asm_acc = ""
        asm_source = ""
        for gid in gids:
            candidate = gene_to_asm.get(gid, "")
            if candidate:
                asm_acc = candidate
                asm_source = f"elink:gene({gid})"
                break

        # Try scaffold fallback
        if not asm_acc:
            scaffold = tx_to_scaffold.get(acc, "")
            if scaffold and scaffold in scaffold_to_asm:
                asm_acc = scaffold_to_asm[scaffold]
                asm_source = f"nuccore:scaffold({scaffold})"

        status = "resolved" if asm_acc else "unresolved"
        rows.append(
            {
                "transcript_id": acc,
                "gene_ids": ";".join(gids),
                "genomic_accession": row.get("genomic_accession", ""),
                "assembly_accession": asm_acc,
                "resolution_source": asm_source,
                "status": status,
                "elink_uids": ";".join(
                    [uid for gid in gids for uid in gene_asm_link.get(gid, [])]
                ),
            }
        )

    df_out = pd.DataFrame(rows)
    resolved = df_out[df_out["status"] == "resolved"]
    unresolved = df_out[df_out["status"] == "unresolved"]

    n_elink = df_out["resolution_source"].str.startswith("elink", na=False).sum()
    n_nuccore = df_out["resolution_source"].str.startswith("nuccore", na=False).sum()

    # Stats relative to unresolved subset only
    unresolved_no_gid = [
        acc for acc in unresolved["transcript_id"] if not tx_to_gids.get(acc)
    ]
    unresolved_no_gid_no_scaffold = [
        acc for acc in unresolved_no_gid if acc not in tx_to_scaffold
    ]

    print(f"\n{'='*55}")
    print(f"Input transcripts             : {len(no_asm)}")
    n_with_gid = len(no_asm) - len(no_geneid_txs)
    print(f"  w/ Gene ID                  : {n_with_gid}")
    print(
        f"  w/ scaffold only            : {len([a for a in no_geneid_txs if a in tx_to_scaffold])}"
    )
    print(
        f"  no Gene ID & no scaffold    : {len([a for a in no_geneid_txs if a not in tx_to_scaffold])}"
    )
    print(f"Resolved via elink            : {n_elink}")
    print(f"Resolved via scaffold fallback: {n_nuccore}")
    print(f"Total resolved                : {len(resolved)}")
    print(f"Still unresolved              : {len(unresolved)}")
    if not unresolved.empty:
        print(f"  of which: no Gene ID        : {len(unresolved_no_gid)}")
        print(
            f"  of which: no GeneID & no scaffold: {len(unresolved_no_gid_no_scaffold)}"
        )

    df_out.to_csv(args.output, sep="\t", index=False)
    print(f"\nWritten → {args.output}")

    # Quick sample of unresolved to diagnose
    if not unresolved.empty:
        print(f"\nFirst 10 unresolved rows:")
        print(
            unresolved[["transcript_id", "gene_ids", "genomic_accession"]]
            .head(10)
            .to_string(index=False)
        )


if __name__ == "__main__":
    main()
