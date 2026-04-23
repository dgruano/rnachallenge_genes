#!/usr/bin/env python3
"""
Extract NW_ genome accessions from local NCBI Gene XML files.
Useful for testing or batch processing downloaded XML responses.
"""

import sys
import os
import xml.etree.ElementTree as ET

def extract_nw_from_file(xml_file):
    """Extract NW_ accession, GeneID, and organism from an XML file."""
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
    except ET.ParseError as e:
        return None, f"Parse error: {e}"

    # Handle Entrezgene-Set wrapper
    if root.tag == "Entrezgene-Set":
        entrezgene = root.find("Entrezgene")
        if entrezgene is None:
            return None, "No Entrezgene element found"
    else:
        entrezgene = root

    result = {}

    # Extract GeneID
    geneid_el = entrezgene.find(".//Gene-track_geneid")
    if geneid_el is not None and geneid_el.text:
        result["gene_id"] = geneid_el.text

    # Extract organism
    org_el = entrezgene.find(".//Org-ref_taxname")
    if org_el is not None and org_el.text:
        result["organism"] = org_el.text

    # Extract NW_ accession (scaffold accession)
    locus_el = entrezgene.find("Entrezgene_locus")
    if locus_el is not None:
        for gc_genomic in locus_el.findall("Gene-commentary"):
            gc_type = gc_genomic.find("Gene-commentary_type")
            if gc_type is not None and gc_type.get("value") == "genomic":
                acc_el = gc_genomic.find("Gene-commentary_accession")
                if acc_el is not None and acc_el.text:
                    result["nw_accession"] = acc_el.text
                    break

    if not result.get("nw_accession"):
        return None, "No NW_ accession found"

    return result, None


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python extract_nw_from_xml.py <xml_file> [xml_file2 ...]")
        sys.exit(1)

    print("GeneID\tOrganism\tNW_ Accession")
    for xml_file in sys.argv[1:]:
        if not os.path.exists(xml_file):
            print(f"ERROR: {xml_file} not found", file=sys.stderr)
            continue

        result, error = extract_nw_from_file(xml_file)
        if error:
            print(f"ERROR in {xml_file}: {error}", file=sys.stderr)
        else:
            print(f"{result.get('gene_id', 'N/A')}\t{result.get('organism', 'N/A')}\t{result.get('nw_accession', 'N/A')}")
