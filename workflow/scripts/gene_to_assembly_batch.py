#!/usr/bin/env python3
"""
Gene ID → Assembly Accession + XM Transcript Resolver  (Batched edition)

Uses NCBI E-utilities with epost/webenv history to minimise request count:

  Phase 1 – Gene info  : epost(gene) + efetch(Entrezgene-Set) per chunk
                         → parses all records in one HTTP response
  Phase 2 – Elink      : elink(gene→assembly) with all IDs comma-separated
                         per chunk → one linkset per gene ID
  Phase 3 – Assembly   : single esummary for all unique assembly UIDs
  Phase 4 – Fallback   : per-gene nuccore→assembly for genes with no elink
                         result (scaffold DBLINK strategy)

Request count (N genes, chunk=200):
  gene_to_assembly.py       ≥ 3 N
  gene_to_assembly_batch.py ≈ 3 * ceil(N/200) + fallback_genes * 3

Requirements:
    pip install biopython
Usage:
    python gene_to_assembly_batch.py --ids 26516672 5715093 --email your@email.com
    python gene_to_assembly_batch.py --file gene_ids.txt --email your@email.com \\
        --api_key YOUR_KEY --output results.tsv
"""

import argparse
import csv
import sys

try:
    from Bio import Entrez
except ImportError:
    print("Error: Biopython is required. Install with: pip install biopython")
    sys.exit(1)

import os as _os

sys.path.insert(0, _os.path.dirname(__file__))
from ncbi_entrez_utils import (
    CHUNK_SIZE,
    batch_fetch_gene_info,
    batch_link_genes_to_assemblies,
    chunks,
    fetch_assembly_from_nuccore,
    resolve_assembly_uids_map,
)

RATE_LIMIT_DELAY = 0.34  # ~3 req/sec without API key
RATE_LIMIT_DELAY_KEY = 0.11  # ~10 req/sec with API key


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args():
    parser = argparse.ArgumentParser(
        description="Retrieve assembly accessions and XM transcripts from NCBI Gene IDs "
        "(batched using epost history)."
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--ids",
        nargs="+",
        metavar="GENE_ID",
        help="One or more NCBI Gene IDs (space-separated)",
    )
    input_group.add_argument(
        "--file", metavar="FILE", help="Text file with one Gene ID per line"
    )
    parser.add_argument("--email", required=True, help="Your email (required by NCBI)")
    parser.add_argument("--api_key", default=None, help="NCBI API key (optional)")
    parser.add_argument(
        "--output", default=None, help="Output TSV file (default: stdout)"
    )
    return parser.parse_args()


def load_gene_ids(args):
    if args.ids:
        return [g.strip() for g in args.ids if g.strip()]
    try:
        with open(args.file) as fh:
            return [l.strip() for l in fh if l.strip() and not l.startswith("#")]
    except FileNotFoundError:
        print(f"Error: File not found: {args.file}", file=sys.stderr)
        sys.exit(1)


# ---------------------------------------------------------------------------
# XML parsing helpers  (moved to ncbi_entrez_utils; kept as reference)
# ---------------------------------------------------------------------------


def _parse_entrezgene_element_local(el):  # kept for any local use; real impl in utils
    """
    Extract structured info from a single parsed <Entrezgene> XML element.

    Returns a dict with keys:
      gene_id, gene_status, locus_tag, scaffold_acc,
      mrna_accs (list), protein_accs (list), update_date
    """
    result = {
        "gene_id": "",
        "gene_status": "unknown",
        "locus_tag": "",
        "scaffold_acc": "",
        "mrna_accs": [],
        "protein_accs": [],
        "update_date": "",
    }

    # Gene ID
    gid_el = el.find(".//Gene-track_geneid")
    if gid_el is not None:
        result["gene_id"] = (gid_el.text or "").strip()

    # Gene status
    status_el = el.find(".//Gene-track_status")
    if status_el is not None:
        val = status_el.get("value", "")
        result["gene_status"] = val if val else status_el.text or "unknown"

    # Locus tag
    lt = el.find(".//Gene-ref_locus-tag")
    if lt is not None:
        result["locus_tag"] = lt.text or ""

    # Update date
    y = el.find(".//Gene-track_update-date//Date-std_year")
    m = el.find(".//Gene-track_update-date//Date-std_month")
    d = el.find(".//Gene-track_update-date//Date-std_day")
    if y is not None:
        result["update_date"] = f"{y.text}/{m.text:>02}/{d.text:>02}"

    # Genomic scaffold + mRNA / protein accessions from Entrezgene_locus
    locus_el = el.find("Entrezgene_locus")
    if locus_el is not None:
        for gc_genomic in locus_el.findall("Gene-commentary"):
            gc_type = gc_genomic.find("Gene-commentary_type")
            if gc_type is None or gc_type.get("value") != "genomic":
                continue

            acc_el = gc_genomic.find("Gene-commentary_accession")
            if acc_el is not None and acc_el.text:
                result["scaffold_acc"] = acc_el.text

            products_el = gc_genomic.find("Gene-commentary_products")
            if products_el is None:
                continue
            for gc_mrna in products_el.findall("Gene-commentary"):
                mrna_type = gc_mrna.find("Gene-commentary_type")
                if mrna_type is None or mrna_type.get("value") != "3":
                    continue
                mrna_acc = gc_mrna.find("Gene-commentary_accession")
                mrna_ver = gc_mrna.find("Gene-commentary_version")
                if mrna_acc is not None and mrna_acc.text:
                    ver = f".{mrna_ver.text}" if mrna_ver is not None else ""
                    result["mrna_accs"].append(f"{mrna_acc.text}{ver}")

                prot_products = gc_mrna.find("Gene-commentary_products")
                if prot_products is None:
                    continue
                for gc_prot in prot_products.findall("Gene-commentary"):
                    prot_type = gc_prot.find("Gene-commentary_type")
                    if prot_type is None or prot_type.get("value") != "8":
                        continue
                    prot_acc = gc_prot.find("Gene-commentary_accession")
                    prot_ver = gc_prot.find("Gene-commentary_version")
                    if prot_acc is not None and prot_acc.text:
                        ver = f".{prot_ver.text}" if prot_ver is not None else ""
                        result["protein_accs"].append(f"{prot_acc.text}{ver}")

    result["mrna_accs"] = list(dict.fromkeys(result["mrna_accs"]))
    result["protein_accs"] = list(dict.fromkeys(result["protein_accs"]))
    return result


# ---------------------------------------------------------------------------
# Batch fetchers  (implementation in ncbi_entrez_utils; imported above)
# ---------------------------------------------------------------------------
#
# batch_fetch_gene_info, batch_link_genes_to_assemblies,
# resolve_assembly_uids_map, fetch_assembly_from_nuccore
# are all imported from ncbi_entrez_utils at the top of this file.

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

FIELDNAMES = [
    "gene_id",
    "gene_status",
    "locus_tag",
    "scaffold_acc",
    "mrna_accessions",
    "protein_accessions",
    "gene_update_date",
    "assembly_accession",
    "assembly_name",
    "seq_release_date",
    "organism",
    "assembly_status",
]

_EMPTY_ASM = {
    "assembly_accession": "",
    "assembly_name": "",
    "seq_release_date": "",
    "organism": "",
    "assembly_status": "",
}


def write_results(rows, output_path):
    if output_path:
        fh = open(output_path, "w", newline="")
    else:
        fh = sys.stdout
    writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    if output_path:
        fh.close()
        print(f"Results written to: {output_path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    args = parse_args()
    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key
    delay = RATE_LIMIT_DELAY_KEY if args.api_key else RATE_LIMIT_DELAY

    gene_ids = load_gene_ids(args)
    if not gene_ids:
        print("Error: No Gene IDs provided.", file=sys.stderr)
        sys.exit(1)

    print(f"Processing {len(gene_ids)} Gene ID(s)...", file=sys.stderr)

    # Phase 1: Batch gene XML fetch
    print("Phase 1: Fetching gene records (batched)...", file=sys.stderr)
    gene_info_map = batch_fetch_gene_info(gene_ids, delay)

    # Phase 2: Batch gene → assembly elink
    print("Phase 2: Linking genes to assemblies (batched)...", file=sys.stderr)
    gene_assembly_map = batch_link_genes_to_assemblies(gene_ids, delay)

    # Phase 3: Batch assembly summary for all unique UIDs from elink
    all_uids = [uid for uids in gene_assembly_map.values() for uid in uids]
    assembly_uid_map = {}
    if all_uids:
        print("Phase 3: Resolving assembly records (batched)...", file=sys.stderr)
        assembly_uid_map = resolve_assembly_uids_map(all_uids, delay)

    # Phase 4: Per-gene fallback + row assembly
    need_fallback = [gid for gid in gene_ids if not gene_assembly_map.get(gid)]
    if need_fallback:
        print(
            f"Phase 4: Nuccore fallback for {len(need_fallback)} gene(s) "
            f"with no assembly link...",
            file=sys.stderr,
        )

    rows = []
    for i, gid in enumerate(gene_ids):
        info = gene_info_map.get(gid, {})
        asm_uids = gene_assembly_map.get(gid, [])

        row = {
            "gene_id": gid,
            "gene_status": info.get("gene_status", ""),
            "locus_tag": info.get("locus_tag", ""),
            "scaffold_acc": info.get("scaffold_acc", ""),
            "mrna_accessions": ";".join(info.get("mrna_accs", [])),
            "protein_accessions": ";".join(info.get("protein_accs", [])),
            "gene_update_date": info.get("update_date", ""),
            **_EMPTY_ASM,
        }

        if asm_uids:
            # Use the first assembly UID returned by elink
            asm = assembly_uid_map.get(asm_uids[0], {})
            row["assembly_accession"] = asm.get("assembly_accession", "N/A")
            row["assembly_name"] = asm.get("assembly_name", "N/A")
            row["seq_release_date"] = asm.get("seq_release_date", "N/A")
            row["organism"] = asm.get("organism", "N/A")
            row["assembly_status"] = asm.get("assembly_status", "N/A")
        elif row["scaffold_acc"]:
            if need_fallback and gid in need_fallback:
                print(
                    f"  [{need_fallback.index(gid)+1}/{len(need_fallback)}] "
                    f"Nuccore fallback for gene {gid}...",
                    file=sys.stderr,
                )
            nuccore_asm = fetch_assembly_from_nuccore(row["scaffold_acc"], delay)
            if nuccore_asm and "error" not in nuccore_asm:
                row["assembly_accession"] = nuccore_asm.get("assembly_accession", "N/A")
                row["assembly_name"] = nuccore_asm.get("assembly_name", "N/A")
                row["seq_release_date"] = nuccore_asm.get("seq_release_date", "N/A")
                row["organism"] = nuccore_asm.get("organism", "N/A")
                row["assembly_status"] = nuccore_asm.get("assembly_status", "N/A")
            else:
                row["assembly_accession"] = "N/A (no assembly link found)"
        else:
            row["assembly_accession"] = "N/A (discontinued, no scaffold)"

        rows.append(row)

    write_results(rows, args.output)
    print(f"Done. {len(rows)} Gene ID(s) resolved.", file=sys.stderr)


if __name__ == "__main__":
    main()
