#!/usr/bin/env python3
"""
Gene ID → Assembly Accession + XM Transcript Resolver
Uses NCBI E-utilities to retrieve:
  - Current assembly linked to each Gene ID (via elink gene→assembly)
  - XM/XR transcript accession(s) from the Gene XML record
  - Genomic scaffold accession (NW_/NC_) from the Gene XML record
  - Gene status (live / discontinued)

NOTE: elink gene→assembly returns the CURRENT (Latest) assembly.
      For discontinued genes the linked assembly may be absent; the
      scaffold accession embedded in the XML (NW_*/NC_*) is the
      historically annotated one and can be used to query assembly_summary
      to identify the original assembly.

Requirements:
    pip install biopython
Usage:
    python gene_to_assembly.py --ids 26516672 5715093 --email your@email.com
    python gene_to_assembly.py --file gene_ids.txt --email your@email.com \
        --api_key YOUR_KEY --output results.tsv
"""

import argparse
import re
import time
import sys
import csv
import xml.etree.ElementTree as ET

try:
    from Bio import Entrez
except ImportError:
    print("Error: Biopython is required. Install with: pip install biopython")
    sys.exit(1)


CHUNK_SIZE = 200
RATE_LIMIT_DELAY     = 0.34   # ~3 req/sec without API key
RATE_LIMIT_DELAY_KEY = 0.11   # ~10 req/sec with API key


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Retrieve assembly accessions and XM transcripts from NCBI Gene IDs."
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--ids",  nargs="+", metavar="GENE_ID",
                              help="One or more NCBI Gene IDs (space-separated)")
    input_group.add_argument("--file", metavar="FILE",
                              help="Text file with one Gene ID per line")
    parser.add_argument("--email",   required=True, help="Your email (required by NCBI)")
    parser.add_argument("--api_key", default=None,  help="NCBI API key (optional)")
    parser.add_argument("--output",  default=None,  help="Output TSV file (default: stdout)")
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


def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


# ---------------------------------------------------------------------------
# XML parsing helpers
# ---------------------------------------------------------------------------

def parse_gene_xml(xml_text):
    """
    Parse a single Entrezgene XML block and extract:
      - gene_status  : 'live' | 'discontinued' | 'unknown'
      - locus_tag    : e.g. CHLREDRAFT_162526
      - scaffold_acc : genomic accession from Entrezgene_locus (e.g. NW_001843471)
      - mrna_accs    : list of XM_/NM_/XR_/NR_ accessions from Entrezgene_locus products
      - protein_accs : list of XP_/NP_ accessions
      - update_date  : last update date string

    XML path for mRNA accession (based on the actual record structure):
      Entrezgene
        └─ Entrezgene_locus
             └─ Gene-commentary  [type genomic]
                  └─ Gene-commentary_products
                       └─ Gene-commentary  [type mRNA, value="3"]
                            └─ Gene-commentary_accession   ← XM accession lives here
    """
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return {}

    ns = {}  # no namespace in these records

    result = {
        "gene_status":   "unknown",
        "locus_tag":     "",
        "scaffold_acc":  "",
        "mrna_accs":     [],
        "protein_accs":  [],
        "update_date":   "",
    }

    # Handle Entrezgene-Set wrapper (e.g., from direct efetch XML file)
    if root.tag == "Entrezgene-Set":
        entrezgene = root.find("Entrezgene")
        if entrezgene is None:
            return result
    else:
        entrezgene = root

    # --- Gene status ---
    status_el = entrezgene.find(".//Gene-track_status")
    if status_el is not None:
        val = status_el.get("value", "")
        result["gene_status"] = val if val else status_el.text or "unknown"

    # --- Locus tag ---
    lt = entrezgene.find(".//Gene-ref_locus-tag")
    if lt is not None:
        result["locus_tag"] = lt.text or ""

    # --- Update date ---
    y = entrezgene.find(".//Gene-track_update-date//Date-std_year")
    m = entrezgene.find(".//Gene-track_update-date//Date-std_month")
    d = entrezgene.find(".//Gene-track_update-date//Date-std_day")
    if y is not None:
        result["update_date"] = f"{y.text}/{m.text:>02}/{d.text:>02}"

    # --- Genomic scaffold + mRNA/protein accessions from Entrezgene_locus ---
    # The structure is:
    # Entrezgene_locus
    #   Gene-commentary (type=1, genomic)        ← scaffold
    #     Gene-commentary_accession              ← NW_ / NC_ accession
    #     Gene-commentary_products
    #       Gene-commentary (type=3, mRNA)       ← XM_ accession
    #         Gene-commentary_accession
    #         Gene-commentary_products
    #           Gene-commentary (type=8, peptide) ← XP_ accession
    #             Gene-commentary_accession

    locus_el = entrezgene.find("Entrezgene_locus")
    if locus_el is not None:
        for gc_genomic in locus_el.findall("Gene-commentary"):
            gc_type = gc_genomic.find("Gene-commentary_type")
            if gc_type is not None and gc_type.get("value") == "genomic":
                # scaffold accession
                acc_el = gc_genomic.find("Gene-commentary_accession")
                if acc_el is not None and acc_el.text:
                    result["scaffold_acc"] = acc_el.text

                # mRNA products of this genomic commentary
                products_el = gc_genomic.find("Gene-commentary_products")
                if products_el is not None:
                    for gc_mrna in products_el.findall("Gene-commentary"):
                        mrna_type = gc_mrna.find("Gene-commentary_type")
                        if mrna_type is not None and mrna_type.get("value") == "3":
                            mrna_acc = gc_mrna.find("Gene-commentary_accession")
                            mrna_ver = gc_mrna.find("Gene-commentary_version")
                            if mrna_acc is not None and mrna_acc.text:
                                ver = f".{mrna_ver.text}" if mrna_ver is not None else ""
                                result["mrna_accs"].append(f"{mrna_acc.text}{ver}")

                            # protein products of this mRNA
                            prot_products = gc_mrna.find("Gene-commentary_products")
                            if prot_products is not None:
                                for gc_prot in prot_products.findall("Gene-commentary"):
                                    prot_type = gc_prot.find("Gene-commentary_type")
                                    if prot_type is not None and prot_type.get("value") == "8":
                                        prot_acc = gc_prot.find("Gene-commentary_accession")
                                        prot_ver = gc_prot.find("Gene-commentary_version")
                                        if prot_acc is not None and prot_acc.text:
                                            ver = f".{prot_ver.text}" if prot_ver is not None else ""
                                            result["protein_accs"].append(f"{prot_acc.text}{ver}")

    result["mrna_accs"]    = list(dict.fromkeys(result["mrna_accs"]))    # deduplicate
    result["protein_accs"] = list(dict.fromkeys(result["protein_accs"]))
    return result


# ---------------------------------------------------------------------------
# Assembly resolution
# ---------------------------------------------------------------------------

# Matches "Assembly: GCF_000002655.1" or "Assembly: GCA_000002655.1" in a
# GenBank flatfile DBLINK block.
_DBLINK_ASSEMBLY_RE = re.compile(r"Assembly:\s*(GC[AF]_\d+\.\d+)", re.IGNORECASE)

# Scaffold accession prefixes that can be searched in nuccore
_SCAFFOLD_PREFIXES = ("NW_", "NC_", "NT_", "NZ_")


def fetch_assembly_accession_from_dblink(nuccore_uid, delay):
    """
    Fetch the GenBank flatfile header for a nuccore UID and parse the
    DBLINK block for an Assembly: accession.  Reading stops at the
    FEATURES line so the (potentially very large) sequence is never
    downloaded.
    Returns an accession string (e.g. 'GCF_000002655.1') or None.
    """
    try:
        handle = Entrez.efetch(db="nuccore", id=nuccore_uid,
                               rettype="gb", retmode="text")
        assembly_acc = None
        for line in handle:
            if line.startswith("FEATURES"):
                break
            m = _DBLINK_ASSEMBLY_RE.search(line)
            if m:
                assembly_acc = m.group(1)
                break
        handle.close()
        time.sleep(delay)
        return assembly_acc
    except Exception:
        return None


def resolve_assembly_by_accession(assembly_acc, delay):
    """
    Resolve a GCA_/GCF_ accession string to full assembly metadata
    via esearch + esummary.
    """
    try:
        search_handle = Entrez.esearch(db="assembly", term=assembly_acc, retmax=1)
        search_result = Entrez.read(search_handle)
        search_handle.close()
        time.sleep(delay)
        uids = search_result.get("IdList", [])
        if not uids:
            return {}
        return resolve_assembly_uids(uids, delay)[0]
    except Exception:
        return {}


def fetch_assembly_from_nuccore(scaffold_acc, delay):
    """
    Use a scaffold accession to query nuccore, then resolve the assembly via:
      1. elink nuccore → assembly  (fast, works for most current records)
      2. DBLINK Assembly: field in the GenBank flatfile header (fallback for
         records where the elink table is not populated)
    Handles NW_, NC_, NT_, NZ_ prefixes.
    """
    if not scaffold_acc or not scaffold_acc.startswith(_SCAFFOLD_PREFIXES):
        return {}

    try:
        # Search nuccore for the scaffold accession
        search_handle = Entrez.esearch(db="nuccore", term=scaffold_acc, retmax=1)
        search_result = Entrez.read(search_handle)
        search_handle.close()
        time.sleep(delay)

        nuccore_uids = search_result.get("IdList", [])
        if not nuccore_uids:
            return {}

        nuccore_uid = nuccore_uids[0]

        # Link nuccore → assembly
        link_handle = Entrez.elink(dbfrom="nuccore", db="assembly", id=nuccore_uid)
        link_result = Entrez.read(link_handle)
        link_handle.close()
        time.sleep(delay)

        assembly_uids = []
        for linkset in link_result:
            for lsd in linkset.get("LinkSetDb", []):
                assembly_uids += [l["Id"] for l in lsd.get("Link", [])]

        if assembly_uids:
            asm_records = resolve_assembly_uids(assembly_uids, delay)
            return asm_records[0] if asm_records else {}

        # Strategy 2: parse DBLINK Assembly: from the GenBank flatfile header
        dblink_acc = fetch_assembly_accession_from_dblink(nuccore_uid, delay)
        if dblink_acc:
            asm = resolve_assembly_by_accession(dblink_acc, delay)
            if asm:
                return asm

        return {}
    except Exception as e:
        return {"error": f"nuccore lookup failed: {e}"}


def resolve_assembly_uids(assembly_uids, delay):
    records = []
    try:
        summary_handle = Entrez.esummary(
            db="assembly", id=",".join(set(assembly_uids))
        )
        summary_result = Entrez.read(summary_handle, validate=False)
        summary_handle.close()
        time.sleep(delay)

        doc_summaries = (summary_result
                         .get("DocumentSummarySet", {})
                         .get("DocumentSummary", []))
        for doc in doc_summaries:
            records.append({
                "assembly_accession": doc.get("AssemblyAccession", "N/A"),
                "assembly_name":      doc.get("AssemblyName",      "N/A"),
                "seq_release_date":   doc.get("SeqReleaseDate",    "N/A"),
                "organism":           doc.get("Organism",          "N/A"),
                "assembly_status":    doc.get("AssemblyStatus",    "N/A"),
            })
    except Exception as e:
        records.append({"assembly_accession": f"ERROR: {e}",
                        "assembly_name": "", "seq_release_date": "",
                        "organism": "", "assembly_status": ""})
    return records or [{"assembly_accession": "N/A", "assembly_name": "N/A",
                         "seq_release_date": "N/A", "organism": "N/A",
                         "assembly_status": "N/A"}]


def fetch_data_for_gene_id(gid, delay):
    """Fetch Gene XML + linked assembly for a single Gene ID."""
    row = {
        "gene_id":           gid,
        "gene_status":       "",
        "locus_tag":         "",
        "scaffold_acc":      "",
        "mrna_accessions":   "",
        "protein_accessions":"",
        "gene_update_date":  "",
        "assembly_accession":"",
        "assembly_name":     "",
        "seq_release_date":  "",
        "organism":          "",
        "assembly_status":   "",
    }

    # 1. Fetch Gene XML to extract transcript accessions
    try:
        fetch_handle = Entrez.efetch(db="gene", id=gid, rettype="xml", retmode="xml")
        xml_text = fetch_handle.read()
        fetch_handle.close()
        time.sleep(delay)

        gene_info = parse_gene_xml(xml_text)
        row["gene_status"]        = gene_info.get("gene_status", "")
        row["locus_tag"]          = gene_info.get("locus_tag", "")
        row["scaffold_acc"]       = gene_info.get("scaffold_acc", "")
        row["mrna_accessions"]    = ";".join(gene_info.get("mrna_accs", []))
        row["protein_accessions"] = ";".join(gene_info.get("protein_accs", []))
        row["gene_update_date"]   = gene_info.get("update_date", "")
    except Exception as e:
        row["mrna_accessions"] = f"XML_ERROR: {e}"

    # 2. Link Gene ID → assembly (via elink or nuccore fallback)
    try:
        link_handle = Entrez.elink(dbfrom="gene", db="assembly", id=gid)
        link_result = Entrez.read(link_handle)
        link_handle.close()
        time.sleep(delay)

        assembly_uids = []
        for linkset in link_result:
            for lsd in linkset.get("LinkSetDb", []):
                assembly_uids += [l["Id"] for l in lsd["Link"]]

        if assembly_uids:
            asm_records = resolve_assembly_uids(assembly_uids, delay)
            # Keep the first (there is usually one reference assembly per gene)
            asm = asm_records[0]
            row["assembly_accession"] = asm["assembly_accession"]
            row["assembly_name"]      = asm["assembly_name"]
            row["seq_release_date"]   = asm["seq_release_date"]
            row["organism"]           = asm["organism"]
            row["assembly_status"]    = asm["assembly_status"]
        else:
            # Fallback: if no elink result but we have scaffold accession,
            # use nuccore → assembly to find the historical assembly
            if row["scaffold_acc"]:
                nuccore_asm = fetch_assembly_from_nuccore(row["scaffold_acc"], delay)
                if nuccore_asm and "error" not in nuccore_asm:
                    row["assembly_accession"] = nuccore_asm.get("assembly_accession", "N/A")
                    row["assembly_name"]      = nuccore_asm.get("assembly_name", "N/A")
                    row["seq_release_date"]   = nuccore_asm.get("seq_release_date", "N/A")
                    row["organism"]           = nuccore_asm.get("organism", "N/A")
                    row["assembly_status"]    = nuccore_asm.get("assembly_status", "N/A")
                else:
                    row["assembly_accession"] = "N/A (no assembly link found)"
            else:
                row["assembly_accession"] = "N/A (discontinued, no scaffold)"
    except Exception as e:
        row["assembly_accession"] = f"ELINK_ERROR: {e}"

    return row


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

FIELDNAMES = [
    "gene_id", "gene_status", "locus_tag",
    "scaffold_acc", "mrna_accessions", "protein_accessions",
    "gene_update_date",
    "assembly_accession", "assembly_name", "seq_release_date",
    "organism", "assembly_status",
]


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

    rows = []
    for i, gid in enumerate(gene_ids):
        print(f"  [{i+1}/{len(gene_ids)}] Gene ID {gid}...", file=sys.stderr)
        rows.append(fetch_data_for_gene_id(gid, delay))

    write_results(rows, args.output)
    print(f"Done. {len(rows)} Gene ID(s) resolved.", file=sys.stderr)


if __name__ == "__main__":
    main()
