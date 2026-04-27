"""
ncbi_assembly_utils.py
Shared utilities for NCBI assembly and GTF operations
=====================================================

Extracted from resolve_abandoned_accessions.py to avoid duplication.
Used by both resolve_abandoned_accessions and resolve_ncbi_assembly_accessions.

Functions:
  - resolve_assembly_ftp: Map assembly accessions → FTP paths
  - download_gtf: Download and cache GTF files
  - extract_all_from_gtf: Extract coordinates by transcript ID
  - extract_annotations_by_geneid: Extract coordinates by NCBI Gene ID
  - map_genomic_to_assembly_elink: Map genomic accessions (NC_/NT_/NW_) → assembly GCF_/GCA_
"""

import gzip
import logging
import re
import time
import urllib.request
from pathlib import Path
from typing import Any, Optional

from Bio import Entrez
import pandas as pd

# ── Constants ────────────────────────────────────────────────────────────────
_GENEID_ATTR_RE = re.compile(r'GeneID[="](\d+)', re.IGNORECASE)
_BATCH_SIZE = 50
_RATE_LIMIT_DELAY = 0.02


def set_entrez_credentials(email: str, api_key: Optional[str] = None):
    """Configure NCBI Entrez credentials globally."""
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key


def _parse_gtf_attr(attribute: str, key: str) -> str:
    """Extract a quoted value from a GTF attribute string."""
    match = re.search(rf'{re.escape(key)} "([^"]+)"', attribute)
    return match.group(1) if match else ""


# ── Batch-efficient NCBI API wrapper ─────────────────────────────────────────

def _retry_ncbi_call(fn, label: str, max_retries: int = 3, retry_wait: float = 0.5):
    """Call fn(); retry on exception with exponential backoff."""
    for attempt in range(1, max_retries + 1):
        try:
            return fn()
        except Exception as exc:
            if attempt < max_retries:
                wait_time = retry_wait * attempt
                time.sleep(wait_time)
            else:
                raise RuntimeError(f"{label} failed after {max_retries} attempts: {exc}")


# ── Assembly accession → FTP path resolution ─────────────────────────────────

def resolve_assembly_ftp(
    assembly_accessions: list[str],
    log: Optional[logging.Logger] = None,
    max_retries: int = 3,
    retry_wait: float = 0.5,
) -> dict[str, dict]:
    """
    Map assembly accessions to {gtf_url, organism}.

    Individual ``esearch`` calls (one per accession) to collect UIDs,
    then one batch ``esummary`` to resolve all FTP paths at once.

    Parameters
    ----------
    assembly_accessions : list of str
        NCBI assembly accessions (GCF_/GCA_)
    log : logging.Logger, optional
        Logger for informational output
    max_retries : int
        Number of retries for API calls
    retry_wait : float
        Initial wait time between retries (exponential backoff)

    Returns
    -------
    dict[str, dict]
        {assembly_accession: {urls: [list], organism: str}}
    """
    if log is None:
        log = logging.getLogger(__name__)

    # Phase 1: esearch per assembly to collect UIDs
    asm_to_uid: dict[str, str] = {}
    for asm in assembly_accessions:
        def _search(a=asm):
            handle = Entrez.esearch(
                db="assembly", term=f"{a}[Assembly Accession]", retmax=1
            )
            result = Entrez.read(handle)
            handle.close()
            return result

        try:
            search = _retry_ncbi_call(_search, f"esearch(assembly) {asm}", max_retries, retry_wait)
        except RuntimeError as exc:
            log.error(str(exc))
            continue

        if not search["IdList"]:
            log.warning(f"  Assembly not found in NCBI: {asm}")
            continue
        asm_to_uid[asm] = search["IdList"][0]
        time.sleep(_RATE_LIMIT_DELAY)

    if not asm_to_uid:
        return {}

    # Phase 2: one batched esummary for all UIDs
    uid_list = list(set(asm_to_uid.values()))
    uid_to_doc: dict[str, Any] = {}

    for i in range(0, len(uid_list), _BATCH_SIZE):
        chunk = uid_list[i : i + _BATCH_SIZE]

        def _summary(c=chunk):
            handle = Entrez.esummary(db="assembly", id=",".join(c), report="full")
            summary = Entrez.read(handle)
            handle.close()
            return summary

        try:
            summary = _retry_ncbi_call(
                _summary,
                f"esummary(assembly) chunk {i // _BATCH_SIZE + 1}",
                max_retries,
                retry_wait,
            )
        except RuntimeError as exc:
            log.error(str(exc))
            continue

        for doc in summary["DocumentSummarySet"]["DocumentSummary"]:
            uid = doc.attributes.get("uid", "")
            if uid:
                uid_to_doc[uid] = doc
        time.sleep(_RATE_LIMIT_DELAY)

    # Phase 3: build final map
    results: dict[str, dict] = {}
    for asm, uid in asm_to_uid.items():
        doc = uid_to_doc.get(uid)
        if doc is None:
            log.warning(f"  No esummary doc for assembly UID {uid} ({asm})")
            continue
        ftp_path = doc.get("FtpPath_RefSeq") or doc.get("FtpPath_GenBank") or ""
        organism = doc.get("Organism", "")
        if not ftp_path or ftp_path == "na":
            log.warning(f"  No FTP path for assembly {asm}")
            continue
        prefix = ftp_path.rsplit("/", 1)[-1]
        results[asm] = {
            "urls": [
                f"{ftp_path}/{prefix}_genomic.gtf.gz",
                f"{ftp_path}/{prefix}_genomic.gff.gz",
            ],
            "organism": organism,
        }

    return results


# ── GTF download and caching ────────────────────────────────────────────────

def download_gtf(
    assembly_acc: str,
    urls: list[str],
    cache_dir: Path,
    log: Optional[logging.Logger] = None,
    max_retries: int = 3,
    retry_wait: float = 0.5,
) -> Optional[Path]:
    """
    Download the assembly GTF/GFF to cache; return path to ``.gz`` or ``None``.

    Tries each URL in ``urls`` sequentially. Caches under ``{cache_dir}/{assembly_acc}/genomic.gtf.gz``.
    Skips download if the file already exists.

    Parameters
    ----------
    assembly_acc : str
        Assembly accession (used for directory name)
    urls : list of str
        GTF/GFF URLs to try in order
    cache_dir : Path
        Root cache directory
    log : logging.Logger, optional
        Logger for informational output
    max_retries : int
        Number of download attempts per URL
    retry_wait : float
        Wait time between retries (exponential backoff)

    Returns
    -------
    Path or None
        Path to cached .gz file, or None if all attempts failed
    """
    if log is None:
        log = logging.getLogger(__name__)

    asm_dir = cache_dir / assembly_acc
    asm_dir.mkdir(parents=True, exist_ok=True)
    gz_path = asm_dir / "genomic.gtf.gz"

    if gz_path.exists():
        log.info(f"  GTF already cached: {gz_path}")
        return gz_path

    for url in urls:
        file_format = "GFF" if ".gff" in url else "GTF"
        log.info(f"  Downloading {file_format} for {assembly_acc}: {url}")
        for attempt in range(1, max_retries + 1):
            try:
                urllib.request.urlretrieve(url, gz_path)
                log.info(f"  Downloaded → {gz_path} ({gz_path.stat().st_size / 1e6:.1f} MB)")
                return gz_path
            except Exception as exc:
                log.warning(
                    f"  {file_format} download attempt {attempt}/{max_retries} failed for {assembly_acc}: {exc}"
                )
                if gz_path.exists():
                    gz_path.unlink()
                if attempt < max_retries:
                    time.sleep(retry_wait * attempt)

        log.warning(f"  All {file_format} download attempts failed for {assembly_acc}, trying next format…")

    log.error(f"  All download attempts failed for {assembly_acc} (tried {len(urls)} format(s))")
    return None


# ── GTF coordinate extraction ──────────────────────────────────────────────

def extract_all_from_gtf(
    gtf_gz: Path,
    transcript_ids: set[str],
    log: Optional[logging.Logger] = None,
) -> dict[str, dict]:
    """
    Extract gene-level annotations for **all** ``transcript_ids`` in one pair of
    sequential passes over the GTF file.

    Pass 1 — find ``gene_id`` (and ``gene_symbol``) for each transcript.
    Pass 2 — find the ``gene`` feature for every discovered gene_id to
              recover ``chrom``, ``start``, ``end``, ``strand``.

    Coordinates follow the pipeline convention: ``start`` is 0-based
    (GTF start − 1), ``end`` is 1-based inclusive (GTF end).

    Returns ``{transcript_id: {gene_id, gene_symbol, chrom, start, end, strand}}``.
    """
    if log is None:
        log = logging.getLogger(__name__)

    # ── Pass 1: transcript_id → (gene_id, gene_symbol) ───────────────────────
    tx_to_gene: dict[str, tuple[str, str]] = {}
    remaining_tx = set(transcript_ids)

    try:
        with gzip.open(gtf_gz, "rt") as fh:
            for line in fh:
                if not remaining_tx:
                    break
                if line.startswith("#") or "transcript_id" not in line:
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                attr = parts[8]
                tx_id = _parse_gtf_attr(attr, "transcript_id")
                if tx_id not in remaining_tx:
                    continue
                gene_id = _parse_gtf_attr(attr, "gene_id")
                if not gene_id:
                    continue
                gene_symbol = (
                    _parse_gtf_attr(attr, "gene_name")
                    or _parse_gtf_attr(attr, "gene")
                    or gene_id
                )
                tx_to_gene[tx_id] = (gene_id, gene_symbol)
                remaining_tx.discard(tx_id)
    except Exception as exc:
        log.warning(f"  GTF pass-1 failed for {gtf_gz}: {exc}")
        return {}

    if remaining_tx:
        log.warning(
            f"  {len(remaining_tx)} transcript(s) not found in GTF: "
            f"{sorted(remaining_tx)[:5]}{'…' if len(remaining_tx) > 5 else ''}"
        )

    if not tx_to_gene:
        return {}

    # ── Pass 2: gene_id → genomic coordinates ────────────────────────────────
    gene_ids_needed = {g for g, _ in tx_to_gene.values()}
    gene_info: dict[str, dict] = {}

    try:
        with gzip.open(gtf_gz, "rt") as fh:
            for line in fh:
                if not gene_ids_needed:
                    break
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] != "gene":
                    continue
                attr = parts[8]
                gene_id = _parse_gtf_attr(attr, "gene_id")
                if gene_id not in gene_ids_needed:
                    continue
                gene_info[gene_id] = {
                    "chrom":       parts[0],
                    "start":       int(parts[3]) - 1,  # GTF 1-based → 0-based
                    "end":         int(parts[4]),       # GTF inclusive → half-open
                    "strand":      parts[6],
                    "gene_symbol": (
                        _parse_gtf_attr(attr, "gene_name")
                        or _parse_gtf_attr(attr, "gene")
                        or gene_id
                    ),
                }
                gene_ids_needed.discard(gene_id)
    except Exception as exc:
        log.warning(f"  GTF pass-2 failed for {gtf_gz}: {exc}")

    if gene_ids_needed:
        log.warning(
            f"  {len(gene_ids_needed)} gene feature(s) not found: "
            f"{sorted(gene_ids_needed)[:5]}{'…' if len(gene_ids_needed) > 5 else ''}"
        )

    # ── Combine ──────────────────────────────────────────────────────────────
    result: dict[str, dict] = {}
    for tx_id, (gene_id, gene_sym) in tx_to_gene.items():
        info = gene_info.get(gene_id)
        if info is None:
            continue
        result[tx_id] = {
            "gene_id":     gene_id,
            "gene_symbol": info["gene_symbol"] or gene_sym,
            "chrom":       info["chrom"],
            "start":       info["start"],
            "end":         info["end"],
            "strand":      info["strand"],
        }

    return result


def extract_annotations_by_geneid(
    gtf_gz: Path,
    geneid_to_transcripts: dict[str, list[str]],
    log: Optional[logging.Logger] = None,
) -> dict[str, dict]:
    """
    Fallback GTF/GFF scan: locate ``gene`` features by NCBI numeric Gene ID.

    Used when the original transcript accession is not found directly,
    but the parent gene is still present in the GTF.

    Parameters
    ----------
    gtf_gz : Path
        Path to compressed GTF/GFF file
    geneid_to_transcripts : dict
        Maps {numeric_gene_id_str: [transcript_id, ...]}
    log : logging.Logger, optional
        Logger for informational output

    Returns
    -------
    dict
        {transcript_id: {gene_id, gene_symbol, chrom, start, end, strand}}
        Caller should set is_ambiguous=True for these rows.
    """
    if log is None:
        log = logging.getLogger(__name__)

    needed: set[str] = set(geneid_to_transcripts.keys())
    gene_features: dict[str, dict] = {}

    try:
        with gzip.open(gtf_gz, "rt") as fh:
            for line in fh:
                if not needed:
                    break
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] != "gene":
                    continue
                attr = parts[8]
                match = _GENEID_ATTR_RE.search(attr)
                if not match:
                    continue
                geneid = match.group(1)
                if geneid not in needed:
                    continue
                gene_features[geneid] = {
                    "chrom":       parts[0],
                    "start":       int(parts[3]) - 1,
                    "end":         int(parts[4]),
                    "strand":      parts[6],
                    "gene_symbol": _parse_gtf_attr(attr, "gene_name") or _parse_gtf_attr(attr, "gene") or geneid,
                }
                needed.discard(geneid)
    except Exception as exc:
        log.warning(f"  GTF GeneID scan failed for {gtf_gz}: {exc}")

    if needed:
        log.warning(f"  {len(needed)} GeneID(s) not found in GTF")

    # Expand: one geneid → multiple transcripts
    result: dict[str, dict] = {}
    for geneid, tx_ids in geneid_to_transcripts.items():
        feature = gene_features.get(geneid)
        if feature is None:
            continue
        for tx_id in tx_ids:
            result[tx_id] = {
                "gene_id":     geneid,
                "gene_symbol": feature["gene_symbol"],
                "chrom":       feature["chrom"],
                "start":       feature["start"],
                "end":         feature["end"],
                "strand":      feature["strand"],
            }

    return result


# ── Genomic accession → parent assembly via elink ────────────────────────────

def map_genomic_to_assembly_elink(
    accessions: list[str],
    log: Optional[logging.Logger] = None,
    max_retries: int = 3,
    retry_wait: float = 0.5,
) -> dict[str, Optional[str]]:
    """
    Map genomic accessions (NC_/NT_/NW_) to parent assembly accessions (GCF_/GCA_)
    using three batched NCBI API calls instead of downloading full GenBank records.

    Phase 1 — esearch(nuccore) per unique accession → nuccore UID
    Phase 2 — batch elink(nuccore→assembly) → assembly UIDs
    Phase 3 — batch esummary(assembly) → GCF_/GCA_ accession strings

    Parameters
    ----------
    accessions : list of str
        Genomic accessions, e.g. ["NC_000001.11", "NT_033779.5"]
    log : logging.Logger, optional
    max_retries : int
    retry_wait : float

    Returns
    -------
    dict[str, Optional[str]]
        {accession: GCF_accession} — value is None when mapping failed.
    """
    if not accessions:
        return {}

    if log is None:
        log = logging.getLogger(__name__)

    result: dict[str, Optional[str]] = {acc: None for acc in accessions}

    # ── Phase 1: accession string → nuccore UID (one esearch per accession) ──
    acc_to_uid: dict[str, str] = {}
    for acc in accessions:
        def _search(a=acc):
            handle = Entrez.esearch(db="nuccore", term=f"{a}[Accession]", retmax=1)
            rec = Entrez.read(handle)
            handle.close()
            return rec

        try:
            rec = _retry_ncbi_call(_search, f"esearch(nuccore) {acc}", max_retries, retry_wait)
        except RuntimeError as exc:
            log.warning(str(exc))
            continue

        if rec["IdList"]:
            acc_to_uid[acc] = rec["IdList"][0]
        else:
            log.debug(f"  {acc}: no nuccore record found")
        time.sleep(_RATE_LIMIT_DELAY)

    if not acc_to_uid:
        log.warning("No nuccore UIDs retrieved — all accessions unresolvable")
        return result

    # ── Phase 2: batch elink nuccore → assembly ────────────────────────────
    nuccore_uid_to_asm_uids: dict[str, list[str]] = {}
    all_nuccore_uids = list(acc_to_uid.values())

    for i in range(0, len(all_nuccore_uids), _BATCH_SIZE):
        chunk = all_nuccore_uids[i : i + _BATCH_SIZE]

        def _elink(c=chunk):
            handle = Entrez.elink(dbfrom="nuccore", db="assembly", id=",".join(c))
            links = Entrez.read(handle)
            handle.close()
            return links

        try:
            links = _retry_ncbi_call(
                _elink,
                f"elink(nuccore→assembly) chunk {i // _BATCH_SIZE + 1}",
                max_retries,
                retry_wait,
            )
        except RuntimeError as exc:
            log.error(str(exc))
            continue

        for linkset in links:
            from_uids = linkset.get("IdList", [])
            asm_uids = [
                link["Id"]
                for lsdb in linkset.get("LinkSetDb", [])
                for link in lsdb.get("Link", [])
            ]
            for fid in from_uids:
                nuccore_uid_to_asm_uids.setdefault(fid, []).extend(asm_uids)
        time.sleep(_RATE_LIMIT_DELAY)

    # ── Phase 3: batch esummary(assembly) → GCF_ accession ────────────────
    all_asm_uids = list({uid for uids in nuccore_uid_to_asm_uids.values() for uid in uids})
    asm_uid_to_gcf: dict[str, str] = {}

    for i in range(0, len(all_asm_uids), _BATCH_SIZE):
        chunk = all_asm_uids[i : i + _BATCH_SIZE]

        def _summary(c=chunk):
            handle = Entrez.esummary(db="assembly", id=",".join(c), report="full")
            summary = Entrez.read(handle)
            handle.close()
            return summary

        try:
            summary = _retry_ncbi_call(
                _summary,
                f"esummary(assembly) chunk {i // _BATCH_SIZE + 1}",
                max_retries,
                retry_wait,
            )
        except RuntimeError as exc:
            log.error(str(exc))
            continue

        for doc in summary["DocumentSummarySet"]["DocumentSummary"]:
            uid = doc.attributes.get("uid", "")
            gcf = doc.get("AssemblyAccession", "")
            if uid and gcf:
                asm_uid_to_gcf[uid] = gcf
        time.sleep(_RATE_LIMIT_DELAY)

    # ── Assemble final result ─────────────────────────────────────────────
    for acc, nuccore_uid in acc_to_uid.items():
        asm_uids = nuccore_uid_to_asm_uids.get(nuccore_uid, [])
        for asm_uid in asm_uids:
            gcf = asm_uid_to_gcf.get(asm_uid)
            if gcf:
                result[acc] = gcf
                log.debug(f"  {acc} → {gcf}")
                break
        if result[acc] is None and asm_uids:
            log.debug(f"  {acc}: assembly UID(s) {asm_uids} had no GCF_ accession in summary")

    mapped = sum(1 for v in result.values() if v is not None)
    log.info(f"map_genomic_to_assembly_elink: {mapped}/{len(accessions)} mapped")
    return result


# ── UCSC→GCF mapping (shared by merge_resolved for noncode_v4/2016) ──────────

EXTENDED_UCSC_TO_GCF: dict[str, str] = {
    "TAIR10":   "GCF_000001735.4",   # Arabidopsis thaliana
    "CE10":     "GCF_000002985.6",   # Caenorhabditis elegans (WBcel235)
    "DM6":      "GCF_000001215.4",   # Drosophila melanogaster
    "RN6":      "GCF_000001895.5",   # Rattus norvegicus
    "MONDOM5":  "GCF_000002295.2",   # Monodelphis domesticus
    "PONABE2":  "GCF_000001545.5",   # Pongo abelii
    "GALGAL4":  "GCF_000002315.6",   # Gallus gallus (GRCg6a)
    "ORNANA1":  "GCF_000002275.2",   # Ornithorhynchus anatinus
    "BOSTAU6":  "GCF_000003055.6",   # Bos taurus (UMD 3.1.1)
    "DANRER10": "GCF_000002035.6",   # Danio rerio (GRCz11)
    # noncode_v4 extras not in the original noncode UCSC map
    "DANRER7":  "GCF_000002035.5",   # Danio rerio (GRCz10)
    "DM3":      "GCF_000001215.3",   # Drosophila melanogaster (BDGP5)
    "GALGAL3":  "GCF_000002315.5",   # Gallus gallus (Gallus_gallus-2.1)
}


def apply_ucsc_to_gcf_mapping(df: pd.DataFrame) -> pd.DataFrame:
    """Replace UCSC assembly names with GCF_ accessions in assembly_accession column."""
    if df.empty or "assembly_accession" not in df.columns:
        return df
    result = df.copy()
    result["assembly_accession"] = result["assembly_accession"].apply(
        lambda v: EXTENDED_UCSC_TO_GCF.get(str(v).strip().upper(), v)
        if pd.notna(v) else v
    )
    return result
