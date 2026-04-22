"""
ncbi_entrez_utils.py
Shared NCBI E-utilities helpers for batch gene/assembly resolution.
====================================================================
Caller must set ``Entrez.email`` and (optionally) ``Entrez.api_key``
before calling any function here.

Functions
---------
chunks                          — generic list chunker
_parse_entrezgene_element       — parse one <Entrezgene> XML element
batch_fetch_gene_info           — epost + efetch gene XML, returns {gid: info}
batch_link_genes_to_assemblies  — elink gene→assembly, returns {gid: [uid]}
resolve_assembly_uids_map       — batched esummary, returns {uid: metadata}
fetch_assembly_accession_from_dblink — parse Assembly: from GenBank DBLINK header
resolve_assembly_by_accession   — esearch + esummary for a GCx_ accession string
resolve_assembly_uids           — thin list-oriented wrapper over resolve_assembly_uids_map
fetch_assembly_from_nuccore     — scaffold accession → assembly metadata
"""

import re
import time
import xml.etree.ElementTree as ET

from Bio import Entrez


CHUNK_SIZE = 200

# Matches "Assembly: GCF_000002655.1" in a GenBank flatfile DBLINK block.
_DBLINK_ASSEMBLY_RE = re.compile(r"Assembly:\s*(GC[AF]_\d+\.\d+)", re.IGNORECASE)

# Scaffold accession prefixes that fetch_assembly_from_nuccore handles.
_SCAFFOLD_PREFIXES = ("NW_", "NC_", "NT_", "NZ_")


# ── Helpers ───────────────────────────────────────────────────────────────────

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def _parse_entrezgene_element(el):
    """
    Extract structured info from a single parsed <Entrezgene> XML element.

    Returns a dict with keys:
      gene_id, gene_status, locus_tag, scaffold_acc,
      mrna_accs (list), protein_accs (list), update_date
    """
    result = {
        "gene_id":      "",
        "gene_status":  "unknown",
        "locus_tag":    "",
        "scaffold_acc": "",
        "mrna_accs":    [],
        "protein_accs": [],
        "update_date":  "",
    }

    gid_el = el.find(".//Gene-track_geneid")
    if gid_el is not None:
        result["gene_id"] = (gid_el.text or "").strip()

    status_el = el.find(".//Gene-track_status")
    if status_el is not None:
        val = status_el.get("value", "")
        result["gene_status"] = val if val else status_el.text or "unknown"

    lt = el.find(".//Gene-ref_locus-tag")
    if lt is not None:
        result["locus_tag"] = lt.text or ""

    y = el.find(".//Gene-track_update-date//Date-std_year")
    m = el.find(".//Gene-track_update-date//Date-std_month")
    d = el.find(".//Gene-track_update-date//Date-std_day")
    if y is not None:
        result["update_date"] = f"{y.text}/{m.text:>02}/{d.text:>02}"

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

    result["mrna_accs"]    = list(dict.fromkeys(result["mrna_accs"]))
    result["protein_accs"] = list(dict.fromkeys(result["protein_accs"]))
    return result


# ── Batch fetchers ────────────────────────────────────────────────────────────

def batch_fetch_gene_info(gene_ids, delay):
    """
    Fetch Entrezgene XML for all gene IDs using epost + efetch per chunk.
    Parses the full Entrezgene-Set response and maps each record back to its
    Gene ID (read from Gene-track_geneid inside the XML).

    Returns: {gene_id_str: info_dict}
    """
    results = {}
    for chunk in chunks(gene_ids, CHUNK_SIZE):
        try:
            post_handle = Entrez.epost(db="gene", id=",".join(chunk))
            post_result = Entrez.read(post_handle)
            post_handle.close()
            time.sleep(delay)

            fetch_handle = Entrez.efetch(
                db="gene",
                webenv=post_result["WebEnv"],
                query_key=post_result["QueryKey"],
                rettype="xml", retmode="xml",
            )
            xml_bytes = fetch_handle.read()
            fetch_handle.close()
            time.sleep(delay)

            try:
                root = ET.fromstring(xml_bytes)
            except ET.ParseError:
                continue

            elements = (root.findall("Entrezgene")
                        if root.tag == "Entrezgene-Set"
                        else [root])
            for el in elements:
                info = _parse_entrezgene_element(el)
                gid = info.get("gene_id", "")
                if gid:
                    results[gid] = info
        except Exception as e:
            for gid in chunk:
                results.setdefault(gid, {"gene_id": gid, "error": str(e)})

    return results


def batch_link_genes_to_assemblies(gene_ids, delay):
    """
    Send all Gene IDs in each chunk as a single comma-separated elink call
    (no epost).  NCBI returns one LinkSet per input ID when IDs are submitted
    this way, preserving the gene→assembly mapping.

    Returns: {gene_id_str: [assembly_uid_str, ...]}
    """
    results = {}
    for chunk in chunks(gene_ids, CHUNK_SIZE):
        try:
            link_handle = Entrez.elink(
                dbfrom="gene", db="assembly", id=",".join(chunk)
            )
            link_result = Entrez.read(link_handle)
            link_handle.close()
            time.sleep(delay)

            for linkset in link_result:
                for gid in linkset.get("IdList", []):
                    uids = [
                        lnk["Id"]
                        for lsd in linkset.get("LinkSetDb", [])
                        for lnk in lsd.get("Link", [])
                    ]
                    results[gid] = uids
        except Exception:
            for gid in chunk:
                results.setdefault(gid, [])

    return results


def resolve_assembly_uids_map(assembly_uids, delay):
    """
    Fetch assembly summaries for all unique UIDs via batched esummary calls.

    Returns: {uid_str: record_dict}
    record_dict keys: assembly_accession, assembly_name, seq_release_date,
                      organism, assembly_status
    """
    unique = list(dict.fromkeys(assembly_uids))   # deduplicate, preserve order
    if not unique:
        return {}

    uid_map = {}
    for chunk in chunks(unique, CHUNK_SIZE):
        try:
            summary_handle = Entrez.esummary(db="assembly", id=",".join(chunk))
            summary_result = Entrez.read(summary_handle, validate=False)
            summary_handle.close()
            time.sleep(delay)

            doc_summaries = (summary_result
                             .get("DocumentSummarySet", {})
                             .get("DocumentSummary", []))
            for doc in doc_summaries:
                uid = doc.attributes.get("uid", "")
                uid_map[uid] = {
                    "assembly_accession": doc.get("AssemblyAccession", "N/A"),
                    "assembly_name":      doc.get("AssemblyName",      "N/A"),
                    "seq_release_date":   doc.get("SeqReleaseDate",    "N/A"),
                    "organism":           doc.get("Organism",          "N/A"),
                    "assembly_status":    doc.get("AssemblyStatus",    "N/A"),
                }
        except Exception:
            pass

    return uid_map


# ── Per-gene nuccore fallback ─────────────────────────────────────────────────

def fetch_assembly_accession_from_dblink(nuccore_uid, delay):
    """
    Fetch the GenBank flatfile header for a nuccore UID and parse the DBLINK
    block for an ``Assembly: GCx_`` accession.  Stops at FEATURES so the
    full sequence is never downloaded.

    Returns an accession string (e.g. ``'GCF_000002655.1'``) or ``None``.
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
        uid_map = resolve_assembly_uids_map(uids, delay)
        return uid_map.get(uids[0], {})
    except Exception:
        return {}


def resolve_assembly_uids(assembly_uids, delay):
    """Single-chunk variant: list of UIDs → list of record dicts."""
    uid_map = resolve_assembly_uids_map(assembly_uids, delay)
    return [uid_map[u] for u in assembly_uids if u in uid_map]


def fetch_assembly_from_nuccore(scaffold_acc, delay):
    """
    Use a scaffold accession to resolve an assembly via:
      1. elink nuccore → assembly  (fast, works for most current records)
      2. DBLINK ``Assembly:`` field in the GenBank flatfile header (fallback)
    Handles NW_, NC_, NT_, NZ_ prefixes.

    Returns a metadata dict (keys: assembly_accession, assembly_name, …)
    or ``{}`` on failure, or ``{"error": "…"}`` on exception.
    """
    if not scaffold_acc or not scaffold_acc.startswith(_SCAFFOLD_PREFIXES):
        return {}

    try:
        search_handle = Entrez.esearch(db="nuccore", term=scaffold_acc, retmax=1)
        search_result = Entrez.read(search_handle)
        search_handle.close()
        time.sleep(delay)

        nuccore_uids = search_result.get("IdList", [])
        if not nuccore_uids:
            return {}

        nuccore_uid = nuccore_uids[0]

        # Strategy 1: elink nuccore → assembly
        link_handle = Entrez.elink(dbfrom="nuccore", db="assembly", id=nuccore_uid)
        link_result = Entrez.read(link_handle)
        link_handle.close()
        time.sleep(delay)

        assembly_uids = [
            lnk["Id"]
            for linkset in link_result
            for lsd in linkset.get("LinkSetDb", [])
            for lnk in lsd.get("Link", [])
        ]
        if assembly_uids:
            records = resolve_assembly_uids(assembly_uids, delay)
            return records[0] if records else {}

        # Strategy 2: parse DBLINK Assembly: from GenBank flatfile header
        dblink_acc = fetch_assembly_accession_from_dblink(nuccore_uid, delay)
        if dblink_acc:
            asm = resolve_assembly_by_accession(dblink_acc, delay)
            if asm:
                return asm

        return {}
    except Exception as e:
        return {"error": f"nuccore lookup failed: {e}"}
