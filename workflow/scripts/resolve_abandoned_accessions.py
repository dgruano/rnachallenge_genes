"""
scripts/resolve_abandoned_accessions.py
Third-pass resolver: abandoned (withdrawn/suppressed) NCBI accessions
======================================================================
Strategy
--------
Follows the notebook pipeline exactly, then applies a Gene-ID-based fallback:

1. ``efetch(db="nucleotide", rettype="gb")`` per transcript to get its
   GenBank record, then extract the parent genomic accession (NC_/NT_/NW_/AC_)
   from ``dbxrefs`` or the ``comment`` annotation field.  Also records NCBI
   Gene IDs (``/db_xref="GeneID:…"``) for use in Step 6.
2. ``efetch(db="nucleotide", rettype="gb")`` per unique genomic accession to
   get the chromosome/scaffold record, then extract the assembly accession
   from ``dbxrefs`` (``Assembly:<acc>`` entries) and Gene IDs.
3. For each unique assembly accession: ``esearch`` to get the assembly UID,
   collect all UIDs, then ONE batch ``esummary(assembly)`` call to recover
   FTP paths and organism names for all assemblies at once.
4. Download assembly GTF (``*_genomic.gtf.gz``) for each unique assembly;
   cache under ``resources/cache/<assembly_accession>/genomic.gtf.gz``.
5. For each assembly, scan the GTF **once** (two passes total) to extract
   annotations for **all** transcripts that belong to that assembly:
     Pass 1 — collect ``(transcript_id → gene_id)`` for the full set.
     Pass 2 — collect ``(gene_id → chrom/start/end/strand)`` for all
               discovered gene_ids.
   This replaces ``2 × N`` per-transcript scans with ``2 × A`` scans
   (A = number of unique assemblies, A ≪ N).
6. **Gene ID fallback** for transcripts still unresolved after Step 5:
     6a — transcript not found in GTF, but assembly already downloaded:
          scan the existing GTF by NCBI Gene ID (``db_xref "GeneID:…"``)
          rather than transcript ID.
     6b — no assembly could be found in Steps 1–2: use Gene IDs to
          identify the current assembly via ``elink gene→assembly``
          (batched, from ``ncbi_entrez_utils``), download the GTF, then
          scan by Gene ID.
   Rows resolved here carry ``is_ambiguous=True`` because coordinates come
   from the gene feature rather than a direct transcript match.
7. Build resolved / unresolved TSVs.

Input
-----
results/ncbi_genbank_unresolved.tsv — only rows with
    reason == withdrawn_or_suppressed

Output schema (resolved)
------------------------
transcript_id | db_source | gene_id | gene_symbol | organism |
assembly_accession | chrom | start | end | strand | is_ambiguous

Output schema (unresolved)
--------------------------
transcript_id | db_source | reason

Usage
-----
Runs under Snakemake OR standalone via argparse:
  snakemake resolve_abandoned_accessions  (Snakemake mode)
  python resolve_abandoned_accessions.py --input <tsv> --output-resolved <tsv> --output-unresolved <tsv> --config <yaml> ...  (CLI mode)
"""

import argparse
import gzip
import json
import re
import sys
import time
import urllib.request
from pathlib import Path
from typing import Any

from Bio import SeqIO

import pandas as pd
import yaml
from Bio import Entrez

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger
from ncbi_entrez_utils import (
    batch_fetch_gene_info,
    batch_link_genes_to_assemblies,
    resolve_assembly_uids_map,
    fetch_assembly_from_nuccore,
)

# ── Auto-detect mode: Snakemake vs CLI ────────────────────────────────────────
_is_snakemake = "snakemake" in dir()

if not _is_snakemake:
    parser = argparse.ArgumentParser(
        description="Resolve abandoned (withdrawn/suppressed) NCBI accessions via GTF lookup"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to ncbi_genbank_unresolved.tsv input file",
    )
    parser.add_argument(
        "--output-resolved",
        required=True,
        help="Output path for resolved IDs",
    )
    parser.add_argument(
        "--output-unresolved",
        required=True,
        help="Output path for unresolved IDs",
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to config.yaml file (YAML or JSON)",
    )
    parser.add_argument(
        "--cache-dir",
        default="resources/cache",
        help="Cache directory for downloaded GTFs (default: resources/cache)",
    )
    parser.add_argument(
        "--log",
        default=None,
        help="Log file path (optional)",
    )

    args = parser.parse_args()

    # Create a mock snakemake namespace from args
    class _Input:
        def __init__(self, unresolved):
            self.unresolved = unresolved

    class _Output:
        def __init__(self, resolved, unresolved):
            self.resolved = resolved
            self.unresolved = unresolved

        def __getitem__(self, idx):
            if idx == 0:
                return self.resolved
            raise IndexError(f"Output index {idx} out of range")

    class _Snakemake:
        def __init__(self, input_path, output_resolved, output_unresolved, config_path, cache_dir, log_path):
            self.input = _Input(input_path)
            self.output = _Output(output_resolved, output_unresolved)
            self.log = [log_path] if log_path else ["/dev/null"]

            # Load config from YAML or JSON
            with open(config_path) as f:
                if config_path.endswith(".json"):
                    self.config = json.load(f)
                else:
                    self.config = yaml.safe_load(f)

            # Override cache_dir if provided via CLI
            self.config["cache_dir"] = args.cache_dir

    snakemake = _Snakemake(  # type: ignore[name-defined]
        args.input,
        args.output_resolved,
        args.output_unresolved,
        args.config,
        args.cache_dir,
        args.log,
    )

# ── Snakemake interface ──────────────────────────────────────────────────────
log = get_logger(
    "resolve_abandoned_accessions", snakemake.log[0]  # type: ignore[name-defined]
)
input_unresolved: str = str(snakemake.input.unresolved)  # type: ignore[name-defined]
out_resolved: str = str(snakemake.output.resolved)  # type: ignore[name-defined]
out_unresolved: str = str(snakemake.output.unresolved)  # type: ignore[name-defined]
cfg: dict[str, Any] = snakemake.config  # type: ignore[name-defined]

Entrez.email = cfg["ncbi_email"]
Entrez.api_key = cfg["ncbi_api_key"]

CACHE_DIR = Path(cfg["cache_dir"])
MAX_RETRIES = int(cfg.get("max_retries", 3))
RETRY_WAIT = float(cfg.get("retry_wait_seconds", 0.5))  # Reduced from 5; exponential backoff: 0.5s, 1s, 1.5s
BATCH_SIZE = int(cfg.get("ncbi_batch_size", 50))
EFETCH_BATCH_SIZE = 50  # Batch up to 50 accessions per efetch call
RATE_LIMIT_DELAY = 0.02  # 50 req/s = 0.02s minimum; API key allows 10 req/s base, but batching helps

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
UNRESOLVED_COLS = ["transcript_id", "db_source", "reason"]


# ── Retry helper ─────────────────────────────────────────────────────────────

def _retry(fn, label: str):
    """Call ``fn()``; retry up to MAX_RETRIES times on any exception."""
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            return fn()
        except Exception as exc:
            log.warning(f"{label} attempt {attempt}/{MAX_RETRIES}: {exc}")
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_WAIT * attempt)
    raise RuntimeError(f"{label} failed after {MAX_RETRIES} attempts")


# ── Step 1: efetch transcript → genomic accession ───────────────────────────
# ── Step 2: efetch genomic record → assembly accession ──────────────────────

_GENOMIC_PREFIXES = ("NC_", "NT_", "NW_", "AC_")
_GENOMIC_RE = re.compile(r'\b((?:NC|NT|NW|AC)_\d+(\.\d+)?)\b')  # Version suffix is not mandatory


def _extract_genomic_accession_from_record(gb_record) -> str | None:
    """Extract the parent genomic (chromosome/scaffold) accession.

    Mirrors the notebook ``extract_genomic_accession`` function: first checks
    ``dbxrefs`` for NC_/NT_/NW_/AC_ prefixes, then falls back to a regex
    search over the ``comment`` annotation field.
    """
    for xref in getattr(gb_record, "dbxrefs", []):
        if xref.startswith(_GENOMIC_PREFIXES):
            return xref
    comment = gb_record.annotations.get("comment", "")
    m = _GENOMIC_RE.search(comment)
    return m.group(1) if m else None


def _extract_assembly_from_genomic_record(gb_record) -> str | None:
    """Extract the assembly accession from a genomic GenBank record.

    Mirrors the notebook ``get_assembly`` function: looks for
    ``Assembly:<acc>`` entries in ``dbxrefs``.
    """
    for xref in getattr(gb_record, "dbxrefs", []):
        if xref.startswith("Assembly:"):
            return xref.split(":", 1)[1]
    return None


def _extract_geneids_from_record(gb_record) -> str | None:
    r"""Extract GeneID numbers from /db_xref qualifiers.

    Searches for /db_xref="GeneID:\d+" in the record qualifiers and returns
    the numeric part(s) as a delimited string. Returns None if no GeneIDs found.
    """
    geneids = []
    for feature in getattr(gb_record, "features", []):
        qualifiers = getattr(feature, "qualifiers", {})
        db_xref = qualifiers.get("db_xref", [])
        for xref in db_xref:
            m = re.search(r"GeneID:(\d+)", xref)
            if m:
                geneids.append(m.group(1))
    return ";".join(geneids) if geneids else None


def _extract_comment_accessions_from_record(gb_record) -> str | None:
    """Extract all genomic accessions (NC_/NT_/NW_/AC_) from comments.

    Searches the comment annotation field for all matching accessions and
    returns them as a delimited string. Returns None if none found.
    """
    comment = gb_record.annotations.get("comment", "")
    matches = _GENOMIC_RE.findall(comment)
    if matches:
        accessions = [m[0] for m in matches]
        return ";".join(accessions)
    return None


def fetch_genomic_accessions(accessions: list[str]) -> tuple[dict[str, str], dict[str, dict]]:
    """Step 1: ``{transcript_acc: genomic_accession}``.

    Batches up to EFETCH_BATCH_SIZE accessions per ``efetch(nucleotide, rettype=gb)`` call
    instead of one per accession. The GenBank record for a withdrawn/suppressed transcript
    still carries the parent genomic accession in its ``dbxrefs`` and ``comment`` fields.

    Returns: (results dict, tracking dict with detailed info per transcript)
    """
    results: dict[str, str] = {}
    tracking: dict[str, dict] = {}  # transcript_id → {dbxrefs, comment, genomic_acc, etc}
    total = len(accessions)

    # Process in batches of EFETCH_BATCH_SIZE
    for batch_start in range(0, len(accessions), EFETCH_BATCH_SIZE):
        batch_end = min(batch_start + EFETCH_BATCH_SIZE, len(accessions))
        batch = accessions[batch_start:batch_end]
        batch_ids = ",".join(batch)

        if batch_end % 50 == 0 or batch_end == total:
            log.info(f"  Step 1: {batch_end}/{total} transcripts fetched (batch {batch_start // EFETCH_BATCH_SIZE + 1})")

        def _fetch(ids=batch_ids):
            handle = Entrez.efetch(
                db="nucleotide", id=ids, rettype="gb", retmode="text"
            )
            return SeqIO.parse(handle, "genbank")

        try:
            gb_records = _retry(_fetch, f"efetch(nucleotide) batch [{batch[0]}...{batch[-1]}]")
            for gb in gb_records:
                acc = gb.id
                if acc not in results:  # Avoid duplicates if parsing returns multiple
                    dbxrefs = getattr(gb, "dbxrefs", [])
                    comment = gb.annotations.get("comment", "")
                    organism = gb.annotations.get("organism", "")
                    genomic = _extract_genomic_accession_from_record(gb)

                    geneids_from_tx = _extract_geneids_from_record(gb)
                    tracking[acc] = {
                        "transcript_id": acc,
                        "transcript_dbxrefs": "|".join(dbxrefs),
                        "transcript_organism": organism,
                        "genomic_accession": genomic or "",
                        "transcript_geneids": geneids_from_tx or "",
                    }

                    if genomic:
                        results[acc] = genomic
                    else:
                        log.debug(f"  {acc}: no genomic accession found in GenBank record")
        except RuntimeError as exc:
            log.warning(f"  Batch [{batch[0]}...{batch[-1]}]: efetch failed — {exc}")

        time.sleep(RATE_LIMIT_DELAY)

    return results, tracking


def fetch_assembly_accessions(genomic_accs: list[str]) -> tuple[dict[str, str], dict[str, dict]]:
    """Step 2: ``{genomic_accession: assembly_accession}``.

    Deduplicated: only unique genomic accessions are fetched.  Batches up to
    EFETCH_BATCH_SIZE per call. The assembly accession is read from the
    ``Assembly:`` entry in ``dbxrefs`` of the genomic GenBank record.

    Returns: (results dict, tracking dict with detailed info per genomic accession)
    """
    unique = list(dict.fromkeys(genomic_accs))  # preserve order, deduplicate
    results: dict[str, str] = {}
    tracking: dict[str, dict] = {}  # genomic_acc → {dbxrefs, organism, assembly_acc, etc}
    total = len(unique)

    # Process in batches of EFETCH_BATCH_SIZE
    for batch_start in range(0, len(unique), EFETCH_BATCH_SIZE):
        batch_end = min(batch_start + EFETCH_BATCH_SIZE, len(unique))
        batch = unique[batch_start:batch_end]
        batch_ids = ",".join(batch)

        if batch_end % 20 == 0 or batch_end == total:
            log.info(f"  Step 2: {batch_end}/{total} genomic records fetched (batch {batch_start // EFETCH_BATCH_SIZE + 1})")

        def _fetch(ids=batch_ids):
            handle = Entrez.efetch(
                db="nucleotide", id=ids, rettype="gb", retmode="text"
            )
            return SeqIO.parse(handle, "genbank")

        try:
            gb_records = _retry(_fetch, f"efetch(nucleotide/genomic) batch [{batch[0]}...{batch[-1]}]")
            for gb in gb_records:
                gacc = gb.id
                if gacc not in results:  # Avoid duplicates if parsing returns multiple
                    dbxrefs = getattr(gb, "dbxrefs", [])
                    organism = gb.annotations.get("organism", "")
                    asm = _extract_assembly_from_genomic_record(gb)
                    geneids = _extract_geneids_from_record(gb)
                    comment_accs = _extract_comment_accessions_from_record(gb)

                    tracking[gacc] = {
                        "genomic_accession": gacc,
                        "genomic_dbxrefs": "|".join(dbxrefs),
                        "genomic_organism": organism,
                        "assembly_accession": asm or "",
                        "genomic_geneids": geneids or "",
                    }

                    if asm:
                        results[gacc] = asm
                    else:
                        log.warning(f"  {gacc}: no Assembly: entry in dbxrefs")
        except RuntimeError as exc:
            log.warning(f"  Batch [{batch[0]}...{batch[-1]}]: efetch failed — {exc}")

        time.sleep(RATE_LIMIT_DELAY)

    return results, tracking





# ── Step 3: Assembly accessions → FTP URLs (batched esummary) ────────────────

def resolve_assembly_ftp(assembly_accessions: list[str]) -> dict[str, dict]:
    """
    Map assembly accessions to ``{gtf_url, organism}``.

    Individual ``esearch`` calls are unavoidable (one per accession) but the
    resulting UIDs are then resolved in one batch ``esummary``, cutting the
    total calls from ``2 × A`` to ``A + ceil(A / BATCH_SIZE)``.
    """
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
            search = _retry(_search, f"esearch(assembly) {asm}")
        except RuntimeError as exc:
            log.error(str(exc))
            continue

        if not search["IdList"]:
            log.warning(f"  Assembly not found in NCBI: {asm}")
            continue
        asm_to_uid[asm] = search["IdList"][0]
        time.sleep(RATE_LIMIT_DELAY)

    if not asm_to_uid:
        return {}

    # Phase 2: one batched esummary for all UIDs
    uid_list = list(set(asm_to_uid.values()))
    uid_to_doc: dict[str, Any] = {}

    for i in range(0, len(uid_list), BATCH_SIZE):
        chunk = uid_list[i : i + BATCH_SIZE]

        def _summary(c=chunk):
            handle = Entrez.esummary(db="assembly", id=",".join(c), report="full")
            summary = Entrez.read(handle)
            handle.close()
            return summary

        try:
            summary = _retry(_summary, f"esummary(assembly) FTP chunk {i // BATCH_SIZE + 1}")
        except RuntimeError as exc:
            log.error(str(exc))
            continue

        for doc in summary["DocumentSummarySet"]["DocumentSummary"]:
            uid = doc.attributes.get("uid", "")
            if uid:
                uid_to_doc[uid] = doc
        time.sleep(RATE_LIMIT_DELAY)

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
        # Approach 3 for future: query NCBI API to discover available formats per assembly
        # (esummary includes FileType info) instead of hardcoding GTF/GFF order
        results[asm] = {
            "urls": [
                f"{ftp_path}/{prefix}_genomic.gtf.gz",
                f"{ftp_path}/{prefix}_genomic.gff.gz",
            ],
            "organism": organism,
        }

    return results


# ── Step 4: Download assembly GTF ────────────────────────────────────────────

def download_gtf(assembly_acc: str, urls: list[str]) -> Path | None:
    """
    Download the assembly GTF/GFF to cache; return path to ``.gz`` or ``None``.

    Tries each URL in ``urls`` sequentially. Caches under ``{CACHE_DIR}/{assembly_acc}/genomic.gtf.gz``.
    Skips the download if the file already exists.
    """
    asm_dir = CACHE_DIR / assembly_acc
    asm_dir.mkdir(parents=True, exist_ok=True)
    gz_path = asm_dir / "genomic.gtf.gz"

    if gz_path.exists():
        log.info(f"  GTF already cached: {gz_path}")
        return gz_path

    for url in urls:
        file_format = "GFF" if ".gff" in url else "GTF"
        log.info(f"  Downloading {file_format} for {assembly_acc}: {url}")
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                urllib.request.urlretrieve(url, gz_path)
                log.info(f"  Downloaded → {gz_path} ({gz_path.stat().st_size / 1e6:.1f} MB)")
                return gz_path
            except Exception as exc:
                log.warning(
                    f"  {file_format} download attempt {attempt}/{MAX_RETRIES} failed for {assembly_acc}: {exc}"
                )
                if gz_path.exists():
                    gz_path.unlink()
                if attempt < MAX_RETRIES:
                    time.sleep(RETRY_WAIT * attempt)

        log.warning(f"  All {file_format} download attempts failed for {assembly_acc}, trying next format…")

    log.error(f"  All download attempts failed for {assembly_acc} (tried {len(urls)} format(s))")
    return None


# ── Step 5: Batch GTF extraction (2 passes per assembly, not per transcript) ──

def _parse_gtf_attr(attribute: str, key: str) -> str:
    """Extract a quoted value from a GTF attribute string."""
    match = re.search(rf'{re.escape(key)} "([^"]+)"', attribute)
    return match.group(1) if match else ""


def extract_all_from_gtf(
    gtf_gz: Path, transcript_ids: set[str]
) -> dict[str, dict]:
    """
    Extract gene-level annotations for **all** ``transcript_ids`` in one
    pair of sequential passes over the GTF file.

    Pass 1 — find ``gene_id`` (and ``gene_symbol``) for each transcript.
    Pass 2 — find the ``gene`` feature for every discovered gene_id to
              recover ``chrom``, ``start``, ``end``, ``strand``.

    Coordinates follow the pipeline convention: ``start`` is 0-based
    (GTF start − 1), ``end`` is 1-based inclusive (GTF end).

    Early-exits each pass as soon as all targets are found.

    Returns ``{transcript_id: {gene_id, gene_symbol, chrom, start, end, strand}}``.
    """
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


# ── Step 5b: Fallback GTF scan by NCBI Gene ID (2 passes replaced by 1) ──────

# Matches GeneID in both GTF (db_xref "GeneID:7157") and GFF3 (Dbxref=GeneID:7157)
_GENEID_ATTR_RE = re.compile(r'GeneID[=:"](\d+)', re.IGNORECASE)


def extract_annotations_by_geneid(
    gtf_gz: Path, geneid_to_transcripts: dict[str, list[str]]
) -> dict[str, dict]:
    """
    Fallback GTF/GFF scan: locate ``gene`` features by NCBI numeric Gene ID
    embedded in ``db_xref "GeneID:<n>"`` (GTF) or ``Dbxref=GeneID:<n>``
    (GFF3) attributes.

    Used when the original transcript accession is no longer annotated in the
    GTF (withdrawn / superseded), but the parent gene is still present.

    Parameters
    ----------
    geneid_to_transcripts : {numeric_gene_id_str: [transcript_id, ...]}
        Maps each NCBI numeric Gene ID to the transcript accession(s) that
        should inherit its genomic coordinates.

    Returns
    -------
    {transcript_id: {gene_id, gene_symbol, chrom, start, end, strand}}
    The caller should set ``is_ambiguous=True`` for these rows because no
    direct transcript-level match was made.
    """
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
                m = _GENEID_ATTR_RE.search(attr)
                if not m:
                    continue
                numeric_gid = m.group(1)
                if numeric_gid not in needed:
                    continue

                # Gene symbol: GTF uses quoted values; GFF3 uses key=value
                gene_id_attr = _parse_gtf_attr(attr, "gene_id")
                if not gene_id_attr:
                    # GFF3: try Name= attribute
                    gn = re.search(r"(?:^|;)Name=([^;]+)", attr)
                    gene_id_attr = gn.group(1) if gn else numeric_gid
                gene_symbol = (
                    _parse_gtf_attr(attr, "gene_name")
                    or _parse_gtf_attr(attr, "gene")
                    or gene_id_attr
                )

                gene_features[numeric_gid] = {
                    "gene_id":     gene_id_attr,
                    "gene_symbol": gene_symbol,
                    "chrom":       parts[0],
                    "start":       int(parts[3]) - 1,  # GTF 1-based → 0-based
                    "end":         int(parts[4]),
                    "strand":      parts[6],
                }
                needed.discard(numeric_gid)
    except Exception as exc:
        log.warning(f"  GTF geneid-scan failed for {gtf_gz}: {exc}")
        return {}

    if needed:
        log.warning(
            f"  {len(needed)} gene(s) not found by Gene ID: "
            f"{sorted(needed)[:5]}{'…' if len(needed) > 5 else ''}"
        )

    result: dict[str, dict] = {}
    for numeric_gid, tx_ids in geneid_to_transcripts.items():
        feat = gene_features.get(numeric_gid)
        if feat is None:
            continue
        for tx_id in tx_ids:
            result[tx_id] = dict(feat)

    return result


# ── Main ─────────────────────────────────────────────────────────────────────

log.info("resolve_abandoned_accessions: GTF-based third-pass resolver (batch mode)")

df_in = pd.read_csv(input_unresolved, sep="\t")
df_work = df_in[
    (df_in["db_source"] == "ncbi") & (df_in["reason"] == "withdrawn_or_suppressed")
].copy()

if df_work.empty:
    log.warning("No withdrawn/suppressed NCBI IDs in input — writing empty outputs")
    pd.DataFrame(columns=RESOLVED_COLS).to_csv(out_resolved, sep="\t", index=False)
    pd.DataFrame(columns=UNRESOLVED_COLS).to_csv(out_unresolved, sep="\t", index=False)
    log.info("resolve_abandoned_accessions complete (nothing to do).")
    sys.exit(0)

accessions: list[str] = df_work["transcript_id"].tolist()
log.info(f"Processing {len(accessions)} withdrawn/suppressed accessions")

# Tracking data for debug output
all_tracking: dict[str, dict] = {}

# ── Step 1: efetch each transcript → genomic accession ──────────────────────
log.info("Step 1: efetch each transcript GenBank record → extract genomic accession …")
log.info(f"  (This makes {len(accessions)} individual API calls — one per transcript)")
acc_to_genomic, tx_tracking = fetch_genomic_accessions(accessions)
all_tracking.update(tx_tracking)
log.info(f"  Genomic accessions found: {len(acc_to_genomic)}/{len(accessions)}")

# ── Step 2: efetch each unique genomic record → assembly accession ───────────
unique_genomic = list(dict.fromkeys(acc_to_genomic.values()))
log.info(
    f"Step 2: efetch {len(unique_genomic)} unique genomic records → extract assembly accession …"
)
genomic_to_assembly, genomic_tracking = fetch_assembly_accessions(unique_genomic)
all_tracking.update(genomic_tracking)
log.info(f"  Assembly accessions found: {len(genomic_to_assembly)}/{len(unique_genomic)}")

# Map each transcript → assembly (via genomic)
acc_to_assembly: dict[str, str] = {}
acc_to_organism: dict[str, str] = {}
for acc in accessions:
    gacc = acc_to_genomic.get(acc)
    if gacc:
        asm = genomic_to_assembly.get(gacc)
        if asm:
            acc_to_assembly[acc] = asm

log.info(f"  Transcripts with assembly accession: {len(acc_to_assembly)}/{len(accessions)}")

# ── Step 3: Resolve FTP paths for unique assemblies (batched esummary) ───────
unique_assemblies = list(set(acc_to_assembly.values()))
log.info(
    f"Step 3: Resolving FTP paths for {len(unique_assemblies)} unique assemblies "
    f"(batched esummary) …"
)
ftp_map = resolve_assembly_ftp(unique_assemblies)
log.info(f"  FTP paths resolved: {len(ftp_map)}/{len(unique_assemblies)}")

# ── Step 4: Download GTFs ─────────────────────────────────────────────────────
log.info("Step 4: Downloading assembly GTFs …")
gtf_paths: dict[str, Path] = {}
for asm in unique_assemblies:
    if asm not in ftp_map:
        log.warning(f"  No FTP path for {asm} — skipping download")
        continue
    path = download_gtf(asm, ftp_map[asm]["urls"])
    if path:
        gtf_paths[asm] = path

log.info(f"  GTFs available: {len(gtf_paths)}/{len(unique_assemblies)}")

# ── Step 5: Batch GTF extraction — 2 passes per assembly, not per transcript ─
log.info("Step 5: Extracting gene annotations from GTFs (batch mode) …")

# Group transcripts by their assembly so each GTF is scanned only once
from collections import defaultdict
asm_to_transcripts: dict[str, set[str]] = defaultdict(set)
for acc in accessions:
    asm = acc_to_assembly.get(acc)
    if asm and asm in gtf_paths:
        asm_to_transcripts[asm].add(acc)

# One extract call per assembly
gtf_annotations: dict[str, dict] = {}  # transcript_id → annotation dict
for asm, tx_set in asm_to_transcripts.items():
    gtf_gz = gtf_paths[asm]
    log.info(f"  Scanning {gtf_gz.name} for {len(tx_set)} transcripts …")
    annotations = extract_all_from_gtf(gtf_gz, tx_set)
    gtf_annotations.update(annotations)
    log.info(f"  → {len(annotations)}/{len(tx_set)} transcripts annotated")

# ── Build resolved / unresolved rows ─────────────────────────────────────────
resolved_rows: list[dict] = []
unresolved_rows: list[dict] = []

for acc in accessions:
    asm = acc_to_assembly.get(acc)
    if not asm:
        unresolved_rows.append(
            {"transcript_id": acc, "db_source": "ncbi", "reason": "no_assembly_found"}
        )
        continue

    if asm not in gtf_paths:
        unresolved_rows.append(
            {"transcript_id": acc, "db_source": "ncbi", "reason": "gtf_download_failed"}
        )
        continue

    annotation = gtf_annotations.get(acc)
    if annotation is None:
        unresolved_rows.append(
            {"transcript_id": acc, "db_source": "ncbi", "reason": "not_found_in_gtf"}
        )
        continue

    organism = ftp_map.get(asm, {}).get("organism", "") or acc_to_organism.get(acc, "")
    resolved_rows.append(
        {
            "transcript_id":      acc,
            "db_source":          "ncbi",
            "gene_id":            annotation["gene_id"],
            "gene_symbol":        annotation["gene_symbol"],
            "organism":           organism,
            "assembly_accession": asm,
            "chrom":              annotation["chrom"],
            "start":              annotation["start"],
            "end":                annotation["end"],
            "strand":             annotation["strand"],
            "is_ambiguous":       False,
        }
    )
# ── Step 6: Gene ID fallback for transcripts still unresolved ────────────────
#
# Two sub-cases:
#   6a — transcript not found in GTF but assembly/GTF already downloaded:
#        scan the existing GTF by Gene ID instead of transcript ID.
#   6b — no assembly found at all: use Gene IDs to identify the assembly,
#        download its GTF, then scan by Gene ID.
#
# Gene IDs are collected from tracking data populated in Steps 1 & 2:
#   transcript_geneids  — extracted from the withdrawn transcript's GB record
#   genomic_geneids     — extracted from the parent genomic scaffold record

log.info("Step 6: Gene ID fallback for remaining unresolved transcripts …")


def _collect_geneids_for_tx(acc: str) -> list[str]:
    """Return deduplicated Gene ID list for a transcript from all_tracking."""
    gids: set[str] = set()
    tx_track = all_tracking.get(acc, {})
    for key in ("transcript_geneids",):
        s = tx_track.get(key, "")
        if s:
            gids.update(g.strip() for g in s.split(";") if g.strip())
    gacc = acc_to_genomic.get(acc)
    if gacc:
        gen_track = all_tracking.get(gacc, {})
        s = gen_track.get("genomic_geneids", "")
        if s:
            gids.update(g.strip() for g in s.split(";") if g.strip())
    return list(gids)


# Partition still-unresolved rows by whether a GTF is already available
still_unresolved_with_gtf: list[tuple[str, list[str], str]] = []  # (acc, gids, asm)
still_unresolved_no_asm:   list[tuple[str, list[str]]]      = []  # (acc, gids)

for row in unresolved_rows:
    acc    = row["transcript_id"]
    reason = row["reason"]
    gids   = _collect_geneids_for_tx(acc)
    if not gids:
        continue  # no Gene ID available — cannot use this fallback
    asm = acc_to_assembly.get(acc)
    if reason == "not_found_in_gtf" and asm and asm in gtf_paths:
        still_unresolved_with_gtf.append((acc, gids, asm))
    elif reason in ("no_assembly_found", "gtf_download_failed"):
        still_unresolved_no_asm.append((acc, gids))

log.info(
    f"  Step 6a candidates (GTF available, tx not found): {len(still_unresolved_with_gtf)}"
)
log.info(
    f"  Step 6b candidates (no assembly yet):              {len(still_unresolved_no_asm)}"
)

# ─ Step 6a: Gene ID scan on already-downloaded GTFs ─────────────────────────
fb_annotations_6a: dict[str, dict] = {}

if still_unresolved_with_gtf:
    asm_to_gid_tx_6a: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    for acc, gids, asm in still_unresolved_with_gtf:
        for gid in gids:
            asm_to_gid_tx_6a[asm][gid].append(acc)

    for asm, gid_tx_map in asm_to_gid_tx_6a.items():
        gtf_gz = gtf_paths[asm]
        n_tx = sum(len(v) for v in gid_tx_map.values())
        log.info(
            f"  Step 6a: scanning {gtf_gz.name} by Gene ID "
            f"({len(gid_tx_map)} gene(s), {n_tx} transcript(s)) …"
        )
        ann = extract_annotations_by_geneid(gtf_gz, gid_tx_map)
        fb_annotations_6a.update(ann)
        log.info(f"  → {len(ann)}/{n_tx} transcript(s) annotated")

# ─ Step 6b: Gene IDs → assembly → download GTF → scan ───────────────────────
fb_annotations_6b: dict[str, dict] = {}
tx_to_gid_asm_6b:  dict[str, tuple[str, str]] = {}  # acc → (gid, asm_acc)

if still_unresolved_no_asm:
    all_gene_ids_6b = list(
        {gid for _acc, gids in still_unresolved_no_asm for gid in gids}
    )
    log.info(
        f"  Step 6b: resolving {len(all_gene_ids_6b)} unique Gene ID(s) → assemblies …"
    )

    # Phase A: batch gene info + elink gene → assembly
    gene_info_map_6b  = batch_fetch_gene_info(all_gene_ids_6b, RATE_LIMIT_DELAY)
    gene_asm_link_6b  = batch_link_genes_to_assemblies(all_gene_ids_6b, RATE_LIMIT_DELAY)

    # Phase B: assembly UIDs → accessions
    all_asm_uids_6b = [uid for uids in gene_asm_link_6b.values() for uid in uids]
    asm_uid_map_6b  = resolve_assembly_uids_map(all_asm_uids_6b, RATE_LIMIT_DELAY) if all_asm_uids_6b else {}

    # Phase C: build gene_id → assembly_accession (with nuccore scaffold fallback)
    gene_to_asm_6b: dict[str, str] = {}
    for gid in all_gene_ids_6b:
        asm_uids = gene_asm_link_6b.get(gid, [])
        if asm_uids:
            asm_info = asm_uid_map_6b.get(asm_uids[0], {})
            asm_acc  = asm_info.get("assembly_accession", "")
            if asm_acc and asm_acc != "N/A":
                gene_to_asm_6b[gid] = asm_acc
        if gid not in gene_to_asm_6b:
            scaffold = gene_info_map_6b.get(gid, {}).get("scaffold_acc", "")
            if scaffold:
                asm_data = fetch_assembly_from_nuccore(scaffold, RATE_LIMIT_DELAY)
                asm_acc  = asm_data.get("assembly_accession", "")
                if asm_acc and asm_acc not in ("N/A", ""):
                    gene_to_asm_6b[gid] = asm_acc

    log.info(
        f"  Step 6b: assembly resolved for {len(gene_to_asm_6b)}/{len(all_gene_ids_6b)} gene(s)"
    )

    # Map each transcript to the first Gene ID with a resolved assembly
    for acc, gids in still_unresolved_no_asm:
        for gid in gids:
            asm_acc = gene_to_asm_6b.get(gid)
            if asm_acc:
                tx_to_gid_asm_6b[acc] = (gid, asm_acc)
                break

    # Resolve FTP paths for assemblies not yet in ftp_map
    new_asms_6b = list({asm for _, asm in tx_to_gid_asm_6b.values()} - set(ftp_map.keys()))
    if new_asms_6b:
        log.info(f"  Step 6b: resolving FTP paths for {len(new_asms_6b)} new assembly(ies) …")
        ftp_map.update(resolve_assembly_ftp(new_asms_6b))

    # Download any new GTFs
    for _acc, (gid, asm_acc) in tx_to_gid_asm_6b.items():
        if asm_acc not in gtf_paths and asm_acc in ftp_map:
            path = download_gtf(asm_acc, ftp_map[asm_acc]["urls"])
            if path:
                gtf_paths[asm_acc] = path

    # Scan GTFs by Gene ID (group by assembly to avoid redundant scans)
    asm_to_gid_tx_6b: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    for acc, (gid, asm_acc) in tx_to_gid_asm_6b.items():
        if asm_acc in gtf_paths:
            asm_to_gid_tx_6b[asm_acc][gid].append(acc)

    for asm_acc, gid_tx_map in asm_to_gid_tx_6b.items():
        gtf_gz = gtf_paths[asm_acc]
        n_tx = sum(len(v) for v in gid_tx_map.values())
        log.info(
            f"  Step 6b: scanning {gtf_gz.name} by Gene ID "
            f"({len(gid_tx_map)} gene(s), {n_tx} transcript(s)) …"
        )
        ann = extract_annotations_by_geneid(gtf_gz, gid_tx_map)
        fb_annotations_6b.update(ann)
        log.info(f"  → {len(ann)}/{n_tx} transcript(s) annotated")

# ─ Merge Step 6 fallback results into resolved/unresolved lists ──────────────
fb_annotations = {**fb_annotations_6a, **fb_annotations_6b}

# Build a combined acc → (gid, asm_acc) lookup for Step 6 hits
tx_to_gid_asm_all: dict[str, tuple[str, str]] = {}
for acc, gids, asm in still_unresolved_with_gtf:
    tx_to_gid_asm_all[acc] = (gids[0] if gids else "", asm)
tx_to_gid_asm_all.update(tx_to_gid_asm_6b)

new_unresolved_rows: list[dict] = []
for row in unresolved_rows:
    acc        = row["transcript_id"]
    annotation = fb_annotations.get(acc)
    if annotation is None:
        new_unresolved_rows.append(row)
        continue
    gid_asm  = tx_to_gid_asm_all.get(acc)
    asm_acc  = gid_asm[1] if gid_asm else (acc_to_assembly.get(acc) or "")
    organism = ftp_map.get(asm_acc, {}).get("organism", "")
    resolved_rows.append(
        {
            "transcript_id":      acc,
            "db_source":          "ncbi",
            "gene_id":            annotation["gene_id"],
            "gene_symbol":        annotation["gene_symbol"],
            "organism":           organism,
            "assembly_accession": asm_acc,
            "chrom":              annotation["chrom"],
            "start":              annotation["start"],
            "end":                annotation["end"],
            "strand":             annotation["strand"],
            "is_ambiguous":       True,   # resolved via Gene ID, not transcript
        }
    )

unresolved_rows = new_unresolved_rows
fb_total = len(fb_annotations_6a) + len(fb_annotations_6b)
log.info(f"  Step 6 total: {fb_total} additional transcript(s) resolved via Gene ID fallback")
# ── Write outputs ─────────────────────────────────────────────────────────────
df_unresolved = pd.DataFrame(unresolved_rows, columns=UNRESOLVED_COLS)

df_resolved.to_csv(out_resolved, sep="\t", index=False)
df_unresolved.to_csv(out_unresolved, sep="\t", index=False)

# ── Write detailed debugging TSV with all collected metadata ──────────────────
debug_rows = []
for acc in accessions:
    row = {"transcript_id": acc}

    # Data from transcript record (Step 1)
    if acc in all_tracking:
        row.update({k: v for k, v in all_tracking[acc].items() if k != "transcript_id"})

    # Data from genomic record (Step 2)
    genomic_acc = acc_to_genomic.get(acc)
    if genomic_acc and genomic_acc in all_tracking:
        row.update({k: v for k, v in all_tracking[genomic_acc].items() if k != "genomic_accession"})

    # Assembly and outcome
    assembly_acc = acc_to_assembly.get(acc)
    row["assembly_accession_found"] = assembly_acc or ""

    # Final resolution status
    resolution = "unresolved"
    reason = ""
    if acc in [r["transcript_id"] for r in resolved_rows]:
        resolution = "resolved"
    elif acc in [r["transcript_id"] for r in unresolved_rows]:
        reason = next(r["reason"] for r in unresolved_rows if r["transcript_id"] == acc)
        resolution = "unresolved"

    row["resolution_status"] = resolution
    row["unresolved_reason"] = reason

    debug_rows.append(row)

# Write debug file (output dir + _debug.tsv)
debug_path = out_resolved.replace(".tsv", "_debug.tsv")
if debug_rows:
    df_debug = pd.DataFrame(debug_rows)
    df_debug.to_csv(debug_path, sep="\t", index=False)
    log.info(f"Written → {debug_path}")

log.info("=" * 60)
log.info(f"Input accessions     : {len(accessions)}")
log.info(f"Resolved via GTF     : {len(df_resolved)}")
log.info(f"Still unresolved     : {len(df_unresolved)}")
if not df_unresolved.empty:
    for reason, grp in df_unresolved.groupby("reason"):
        log.info(f"  {reason:<30}: {len(grp)}")
log.info(f"Written → {out_resolved}")
log.info(f"Debug TSV → {debug_path}")
log.info("resolve_abandoned_accessions complete.")
