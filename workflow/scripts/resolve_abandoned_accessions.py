"""
scripts/resolve_abandoned_accessions.py
Third-pass resolver: abandoned (withdrawn/suppressed) NCBI accessions
======================================================================
Strategy — 4-phase pipeline
----------------------------

Phase 1 — Assembly resolution  (all NCBI API calls, no file I/O)

  Step 1. ``efetch(db="nucleotide", rettype="gb")`` per transcript to get its
          GenBank record → extract parent genomic accession (NC_/NT_/NW_/AC_)
          and NCBI Gene IDs (``/db_xref="GeneID:…"``).
  Step 2. ``efetch(db="nucleotide", rettype="gb")`` per unique genomic
          accession → extract assembly accession from dbxrefs and Gene IDs.
  Step 3. For transcripts still lacking an assembly after Steps 1–2, use
          stored Gene IDs:
            A. Batch ``elink gene→assembly`` to get assembly UIDs.
            B. Batch ``esummary(assembly)`` to convert UIDs → accessions.
            C. For genes with no elink hit, try ``fetch_assembly_from_nuccore``
               using the genomic scaffold already retrieved in Step 2.
          No ``batch_fetch_gene_info`` call — scaffold already known from
          tracking; this avoids large XML downloads for discontinued genes.

Phase 2 — GTF acquisition  (unified download pass)

  Step 4. Resolve FTP paths for ALL assemblies (Steps 1–3) in one batched
          ``esummary`` call via ``resolve_assembly_ftp``.
  Step 5. Download all assembly GTFs; cache under
          ``resources/cache/<assembly_accession>/genomic.gtf.gz``.

Phase 3 — GTF search  (no API calls)

  Step 6. For each assembly, scan its GTF twice:
            Strategy A — ``extract_all_from_gtf``: match by transcript ID
                          → ``is_ambiguous=False``.
            Strategy B — ``extract_annotations_by_geneid``: for transcripts
                          not found in Strategy A, scan by Gene ID
                          → ``is_ambiguous=True``.

Phase 4 — Report

  Write resolved TSV, unresolved TSV, and a detailed debug TSV.

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
  python resolve_abandoned_accessions.py --input <tsv> --output-resolved <tsv> \
      --output-unresolved <tsv> --config <yaml> ...  (CLI mode)
"""

import argparse
import gzip
import json
import re
import sys
import time
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from Bio import Entrez, SeqIO

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger
from ncbi_entrez_utils import (
    batch_link_genes_to_assemblies,
    fetch_assembly_from_nuccore,
    resolve_assembly_uids_map,
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
        def __init__(
            self,
            input_path,
            output_resolved,
            output_unresolved,
            config_path,
            cache_dir,
            log_path,
        ):
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
RETRY_WAIT = float(
    cfg.get("retry_wait_seconds", 0.5)
)  # Reduced from 5; exponential backoff: 0.5s, 1s, 1.5s
BATCH_SIZE = int(cfg.get("ncbi_batch_size", 50))
EFETCH_BATCH_SIZE = 50  # Batch up to 50 accessions per efetch call
RATE_LIMIT_DELAY = (
    0.02  # 50 req/s = 0.02s minimum; API key allows 10 req/s base, but batching helps
)

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
_GENOMIC_RE = re.compile(
    r"\b((?:NC|NT|NW|AC)_\d+(\.\d+)?)\b"
)  # Version suffix is not mandatory


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


def fetch_genomic_accessions(
    accessions: list[str],
) -> tuple[dict[str, str], dict[str, dict]]:
    """Step 1: ``{transcript_acc: genomic_accession}``.

    Batches up to EFETCH_BATCH_SIZE accessions per ``efetch(nucleotide, rettype=gb)`` call
    instead of one per accession. The GenBank record for a withdrawn/suppressed transcript
    still carries the parent genomic accession in its ``dbxrefs`` and ``comment`` fields.

    Returns: (results dict, tracking dict with detailed info per transcript)
    """
    results: dict[str, str] = {}
    tracking: dict[str, dict] = (
        {}
    )  # transcript_id → {dbxrefs, comment, genomic_acc, etc}
    total = len(accessions)

    # Process in batches of EFETCH_BATCH_SIZE
    for batch_start in range(0, len(accessions), EFETCH_BATCH_SIZE):
        batch_end = min(batch_start + EFETCH_BATCH_SIZE, len(accessions))
        batch = accessions[batch_start:batch_end]
        batch_ids = ",".join(batch)

        if batch_end % 50 == 0 or batch_end == total:
            log.info(
                f"  Step 1: {batch_end}/{total} transcripts fetched (batch {batch_start // EFETCH_BATCH_SIZE + 1})"
            )

        def _fetch(ids=batch_ids):
            handle = Entrez.efetch(
                db="nucleotide", id=ids, rettype="gb", retmode="text"
            )
            data = handle.read()
            handle.close()
            from io import StringIO

            return list(SeqIO.parse(StringIO(data), "genbank"))

        try:
            gb_records = _retry(
                _fetch, f"efetch(nucleotide) batch [{batch[0]}...{batch[-1]}]"
            )
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
                        log.debug(
                            f"  {acc}: no genomic accession found in GenBank record"
                        )
        except RuntimeError as exc:
            log.warning(f"  Batch [{batch[0]}...{batch[-1]}]: efetch failed — {exc}")

        time.sleep(RATE_LIMIT_DELAY)

    return results, tracking


def fetch_assembly_accessions(
    genomic_accs: list[str],
) -> tuple[dict[str, str], dict[str, dict]]:
    """Step 2: ``{genomic_accession: assembly_accession}``.

    Deduplicated: only unique genomic accessions are fetched.  Batches up to
    EFETCH_BATCH_SIZE per call. The assembly accession is read from the
    ``Assembly:`` entry in ``dbxrefs`` of the genomic GenBank record.

    Returns: (results dict, tracking dict with detailed info per genomic accession)
    """
    unique = list(dict.fromkeys(genomic_accs))  # preserve order, deduplicate
    results: dict[str, str] = {}
    tracking: dict[str, dict] = (
        {}
    )  # genomic_acc → {dbxrefs, organism, assembly_acc, etc}
    total = len(unique)

    # Process in batches of EFETCH_BATCH_SIZE
    for batch_start in range(0, len(unique), EFETCH_BATCH_SIZE):
        batch_end = min(batch_start + EFETCH_BATCH_SIZE, len(unique))
        batch = unique[batch_start:batch_end]
        batch_ids = ",".join(batch)

        if batch_end % EFETCH_BATCH_SIZE == 0 or batch_end == total:
            log.info(
                f"  Step 2: {batch_end}/{total} genomic records fetched (batch {batch_start // EFETCH_BATCH_SIZE + 1})"
            )

        def _fetch(ids=batch_ids):
            handle = Entrez.efetch(
                db="nucleotide", id=ids, rettype="gb", retmode="text"
            )
            data = handle.read()
            handle.close()
            from io import StringIO

            return list(SeqIO.parse(StringIO(data), "genbank"))

        try:
            gb_records = _retry(
                _fetch, f"efetch(nucleotide/genomic) batch [{batch[0]}...{batch[-1]}]"
            )
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
            summary = _retry(
                _summary, f"esummary(assembly) FTP chunk {i // BATCH_SIZE + 1}"
            )
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
                log.info(
                    f"  Downloaded → {gz_path} ({gz_path.stat().st_size / 1e6:.1f} MB)"
                )
                return gz_path
            except Exception as exc:
                log.warning(
                    f"  {file_format} download attempt {attempt}/{MAX_RETRIES} failed for {assembly_acc}: {exc}"
                )
                if gz_path.exists():
                    gz_path.unlink()
                if attempt < MAX_RETRIES:
                    time.sleep(RETRY_WAIT * attempt)

        log.warning(
            f"  All {file_format} download attempts failed for {assembly_acc}, trying next format…"
        )

    log.error(
        f"  All download attempts failed for {assembly_acc} (tried {len(urls)} format(s))"
    )
    return None


# ── Step 5: Batch GTF extraction (2 passes per assembly, not per transcript) ──


def _parse_gtf_attr(attribute: str, key: str) -> str:
    """Extract a quoted value from a GTF attribute string."""
    match = re.search(rf'{re.escape(key)} "([^"]+)"', attribute)
    return match.group(1) if match else ""


def extract_all_from_gtf(gtf_gz: Path, transcript_ids: set[str]) -> dict[str, dict]:
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
                    "chrom": parts[0],
                    "start": int(parts[3]) - 1,  # GTF 1-based → 0-based
                    "end": int(parts[4]),  # GTF inclusive → half-open
                    "strand": parts[6],
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
            "gene_id": gene_id,
            "gene_symbol": info["gene_symbol"] or gene_sym,
            "chrom": info["chrom"],
            "start": info["start"],
            "end": info["end"],
            "strand": info["strand"],
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
                    "gene_id": gene_id_attr,
                    "gene_symbol": gene_symbol,
                    "chrom": parts[0],
                    "start": int(parts[3]) - 1,  # GTF 1-based → 0-based
                    "end": int(parts[4]),
                    "strand": parts[6],
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


# ── Helper: collect Gene IDs for a transcript from tracking dicts ─────────────


def _collect_geneids_for_tx(acc: str) -> list[str]:
    """Return deduplicated Gene ID list from all_tracking for *acc*.

    Prefers Gene IDs extracted from the parent genomic record (more reliable)
    and supplements with those from the transcript record itself.
    """
    gids: set[str] = set()
    s = all_tracking.get(acc, {}).get("transcript_geneids", "")
    if s:
        gids.update(g.strip() for g in s.split(";") if g.strip())
    gacc = acc_to_genomic.get(acc)
    if gacc:
        s = all_tracking.get(gacc, {}).get("genomic_geneids", "")
        if s:
            gids.update(g.strip() for g in s.split(";") if g.strip())
    return list(gids)


# ── Main ─────────────────────────────────────────────────────────────────────

log.info("resolve_abandoned_accessions: GTF-based third-pass resolver (4-phase)")

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

# Tracking data for debug output (populated in Phase 1 Steps 1 & 2)
all_tracking: dict[str, dict] = {}

# ══════════════════════════════════════════════════════════════════════════════
# Phase 1 — Assembly resolution  (all NCBI API calls, no file I/O)
# ══════════════════════════════════════════════════════════════════════════════

# ── Phase 1 / Step 1: transcript → genomic accession ─────────────────────────
log.info("Phase 1 / Step 1: efetch transcript GenBank records → genomic accession …")
acc_to_genomic, tx_tracking = fetch_genomic_accessions(accessions)
all_tracking.update(tx_tracking)
log.info(f"  Genomic accessions found: {len(acc_to_genomic)}/{len(accessions)}")

# ── Phase 1 / Step 2: genomic record → assembly accession ────────────────────
unique_genomic = list(dict.fromkeys(acc_to_genomic.values()))
log.info(
    f"Phase 1 / Step 2: efetch {len(unique_genomic)} unique genomic records "
    f"→ assembly accession …"
)
genomic_to_assembly, genomic_tracking = fetch_assembly_accessions(unique_genomic)
all_tracking.update(genomic_tracking)
log.info(
    f"  Assembly accessions found: {len(genomic_to_assembly)}/{len(unique_genomic)}"
)

# Build transcript → assembly map from Steps 1+2
acc_to_assembly: dict[str, str] = {}
for acc in accessions:
    gacc = acc_to_genomic.get(acc)
    if gacc:
        asm = genomic_to_assembly.get(gacc)
        if asm:
            acc_to_assembly[acc] = asm

log.info(
    f"  Transcripts with assembly accession: {len(acc_to_assembly)}/{len(accessions)}"
)

# ── Phase 1 / Step 3: Gene ID → assembly for transcripts still missing one ───
#
# Collect Gene IDs from tracking (transcript + genomic scaffold records).
# Strategy: batch elink gene→assembly, then batch esummary UIDs→accessions;
# for genes with no elink hit, fall back to fetch_assembly_from_nuccore using
# the genomic scaffold already retrieved in Step 2.
# NOTE: batch_fetch_gene_info is intentionally NOT called here — the scaffold
# is already known from tracking; calling it would trigger large XML downloads
# for discontinued genes and cause hangs.

acc_to_geneids: dict[str, list[str]] = {}
gid_to_accs: dict[str, list[str]] = {}
for acc in accessions:
    if acc not in acc_to_assembly:
        gids = _collect_geneids_for_tx(acc)
        if gids:
            acc_to_geneids[acc] = gids
            for gid in gids:
                gid_to_accs.setdefault(gid, []).append(acc)

all_unresolved_gids = list(gid_to_accs.keys())
n_no_gid = len(accessions) - len(acc_to_assembly) - len(acc_to_geneids)
log.info(
    f"Phase 1 / Step 3: elink gene→assembly for {len(all_unresolved_gids)} unique Gene ID(s) "
    f"({len(acc_to_geneids)} transcripts without assembly, {n_no_gid} with no Gene ID) …"
)

gene_to_asm: dict[str, str] = {}
if all_unresolved_gids:
    # Phase A: batch elink gene → assembly UIDs
    log.info(
        f"  Phase A: elink gene→assembly for {len(all_unresolved_gids)} gene ID(s) …"
    )
    gene_asm_link = batch_link_genes_to_assemblies(
        all_unresolved_gids, RATE_LIMIT_DELAY
    )
    n_linked = sum(1 for uids in gene_asm_link.values() if uids)
    log.info(
        f"  Phase A done: {n_linked}/{len(all_unresolved_gids)} gene(s) linked to ≥1 assembly UID"
    )

    # Phase B: batch esummary UIDs → assembly accessions
    all_asm_uids = [uid for uids in gene_asm_link.values() for uid in uids]
    n_unique_uids = len(set(all_asm_uids))
    log.info(f"  Phase B: esummary for {n_unique_uids} unique assembly UID(s) …")
    asm_uid_map = (
        resolve_assembly_uids_map(all_asm_uids, RATE_LIMIT_DELAY)
        if all_asm_uids
        else {}
    )
    log.info(f"  Phase B done: {len(asm_uid_map)} assembly UID(s) resolved")

    for gid, uids in gene_asm_link.items():
        if uids:
            asm_info = asm_uid_map.get(uids[0], {})
            asm_acc = asm_info.get("assembly_accession", "")
            if asm_acc and asm_acc != "N/A":
                gene_to_asm[gid] = asm_acc

    log.info(
        f"  After Phases A+B: {len(gene_to_asm)}/{len(all_unresolved_gids)} gene(s) → assembly"
    )

    # Phase C: scaffold fallback for genes with no elink result
    gids_no_elink = [gid for gid in all_unresolved_gids if gid not in gene_to_asm]
    if gids_no_elink:
        # Collect ALL unique scaffolds upfront across every unlinked gene, then
        # fetch each scaffold exactly once (success or failure).  The previous
        # approach only cached successes in scaffold_to_asm, so a scaffold that
        # returned empty was retried for every other gene ID mapping to the same
        # transcript — multiplying API calls for discontinued scaffolds.
        _SCAFFOLD_PREFIXES_C = ("NC_", "NT_", "NW_", "NZ_")
        unique_scaffolds: list[str] = list(
            dict.fromkeys(
                scaffold
                for gid in gids_no_elink
                for acc in gid_to_accs.get(gid, [])
                for scaffold in [acc_to_genomic.get(acc, "")]
                if scaffold and scaffold.startswith(_SCAFFOLD_PREFIXES_C)
            )
        )
        n_no_scaffold = len(gids_no_elink) - len(unique_scaffolds)
        log.info(
            f"  Phase C: nuccore fallback for {len(unique_scaffolds)} unique scaffold(s) "
            f"({len(gids_no_elink)} unlinked gene(s); {n_no_scaffold} have no scaffold) …"
        )
        scaffold_to_asm: dict[str, str] = {}
        for i, scaffold in enumerate(unique_scaffolds, 1):
            asm_data = fetch_assembly_from_nuccore(scaffold, RATE_LIMIT_DELAY)
            asm_acc = asm_data.get("assembly_accession", "")
            if asm_acc and asm_acc not in ("N/A", ""):
                scaffold_to_asm[scaffold] = asm_acc
                log.info(f"  [{i}/{len(unique_scaffolds)}] {scaffold} → {asm_acc}")
            else:
                log.info(
                    f"  [{i}/{len(unique_scaffolds)}] {scaffold} → not found"
                    + (f" ({asm_data['error']})" if asm_data.get("error") else "")
                )
        log.info(
            f"  Phase C done: {len(scaffold_to_asm)}/{len(unique_scaffolds)} scaffold(s) resolved"
        )
        for gid in gids_no_elink:
            for acc in gid_to_accs.get(gid, []):
                scaffold = acc_to_genomic.get(acc, "")
                if scaffold and scaffold in scaffold_to_asm:
                    gene_to_asm[gid] = scaffold_to_asm[scaffold]
                    break

    log.info(
        f"  Gene IDs resolved to assembly: {len(gene_to_asm)}/{len(all_unresolved_gids)}"
    )

    # Map each transcript → assembly via its Gene IDs (first match wins)
    for acc, gids in acc_to_geneids.items():
        for gid in gids:
            asm_acc = gene_to_asm.get(gid)
            if asm_acc:
                acc_to_assembly[acc] = asm_acc
                break

log.info(
    f"  Total transcripts with assembly after Step 3: {len(acc_to_assembly)}/{len(accessions)}"
)

# ══════════════════════════════════════════════════════════════════════════════
# Phase 2 — GTF acquisition  (unified download pass — all assemblies at once)
# ══════════════════════════════════════════════════════════════════════════════

unique_assemblies = list(dict.fromkeys(acc_to_assembly.values()))
log.info(
    f"Phase 2 / Step 4: Resolving FTP paths for {len(unique_assemblies)} unique assemblies …"
)
ftp_map = resolve_assembly_ftp(unique_assemblies)
log.info(f"  FTP paths resolved: {len(ftp_map)}/{len(unique_assemblies)}")

log.info("Phase 2 / Step 5: Downloading assembly GTFs …")
gtf_paths: dict[str, Path] = {}
for asm in unique_assemblies:
    if asm not in ftp_map:
        log.warning(f"  No FTP path for {asm} — skipping download")
        continue
    path = download_gtf(asm, ftp_map[asm]["urls"])
    if path:
        gtf_paths[asm] = path

log.info(f"  GTFs available: {len(gtf_paths)}/{len(unique_assemblies)}")

# ══════════════════════════════════════════════════════════════════════════════
# Phase 3 — GTF search  (no API calls)
# ══════════════════════════════════════════════════════════════════════════════

# ── Phase 3 / Step 6: GTF search — Strategy A (transcript ID), then B (Gene ID)
log.info("Phase 3 / Step 6: Scanning GTFs for transcript annotations …")

# Group transcripts by assembly so each GTF is scanned only once
asm_to_transcripts: dict[str, set[str]] = defaultdict(set)
for acc in accessions:
    asm = acc_to_assembly.get(acc)
    if asm and asm in gtf_paths:
        asm_to_transcripts[asm].add(acc)

# Strategy A: scan by transcript ID
gtf_annotations: dict[str, dict] = {}
for asm, tx_set in asm_to_transcripts.items():
    gtf_gz = gtf_paths[asm]
    log.info(f"  [Strategy A] {asm}: scanning for {len(tx_set)} transcript(s) …")
    ann = extract_all_from_gtf(gtf_gz, tx_set)
    gtf_annotations.update(ann)
    log.info(f"  → {len(ann)}/{len(tx_set)} found")

# Strategy B: Gene ID scan for transcripts not found by Strategy A
asm_to_gid_tx: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
for acc in accessions:
    if acc in gtf_annotations:
        continue
    asm = acc_to_assembly.get(acc)
    if not asm or asm not in gtf_paths:
        continue
    for gid in _collect_geneids_for_tx(acc):
        asm_to_gid_tx[asm][gid].append(acc)

gtf_annotations_geneid: dict[str, dict] = {}
for asm, gid_tx_map in asm_to_gid_tx.items():
    gtf_gz = gtf_paths[asm]
    n_tx = sum(len(v) for v in gid_tx_map.values())
    log.info(
        f"  [Strategy B] {asm}: scanning by Gene ID "
        f"({len(gid_tx_map)} gene(s), {n_tx} transcript(s)) …"
    )
    ann = extract_annotations_by_geneid(gtf_gz, gid_tx_map)
    gtf_annotations_geneid.update(ann)
    log.info(f"  → {len(ann)}/{n_tx} found")

# ══════════════════════════════════════════════════════════════════════════════
# Phase 4 — Report
# ══════════════════════════════════════════════════════════════════════════════

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

    organism = ftp_map.get(asm, {}).get("organism", "")

    # Strategy A — transcript-level match
    annotation = gtf_annotations.get(acc)
    if annotation:
        resolved_rows.append(
            {
                "transcript_id": acc,
                "db_source": "ncbi",
                "gene_id": annotation["gene_id"],
                "gene_symbol": annotation["gene_symbol"],
                "organism": organism,
                "assembly_accession": asm,
                "chrom": annotation["chrom"],
                "start": annotation["start"],
                "end": annotation["end"],
                "strand": annotation["strand"],
                "is_ambiguous": False,
            }
        )
        continue

    # Strategy B — Gene ID match
    annotation = gtf_annotations_geneid.get(acc)
    if annotation:
        resolved_rows.append(
            {
                "transcript_id": acc,
                "db_source": "ncbi",
                "gene_id": annotation["gene_id"],
                "gene_symbol": annotation["gene_symbol"],
                "organism": organism,
                "assembly_accession": asm,
                "chrom": annotation["chrom"],
                "start": annotation["start"],
                "end": annotation["end"],
                "strand": annotation["strand"],
                "is_ambiguous": True,  # resolved via Gene ID, not transcript
            }
        )
        continue

    unresolved_rows.append(
        {"transcript_id": acc, "db_source": "ncbi", "reason": "not_found_in_gtf"}
    )

n_strat_a = sum(1 for r in resolved_rows if not r["is_ambiguous"])
n_strat_b = sum(1 for r in resolved_rows if r["is_ambiguous"])
log.info(
    f"  Strategy A (transcript ID): {n_strat_a} resolved; "
    f"Strategy B (Gene ID): {n_strat_b} resolved"
)

# ── Write outputs ─────────────────────────────────────────────────────────────
df_resolved = pd.DataFrame(resolved_rows, columns=RESOLVED_COLS)
df_unresolved = pd.DataFrame(unresolved_rows, columns=UNRESOLVED_COLS)

df_resolved.to_csv(out_resolved, sep="\t", index=False)
df_unresolved.to_csv(out_unresolved, sep="\t", index=False)

# ── Write detailed debugging TSV ──────────────────────────────────────────────
resolved_set = {r["transcript_id"] for r in resolved_rows}
unresolved_map = {r["transcript_id"]: r["reason"] for r in unresolved_rows}

debug_rows = []
for acc in accessions:
    row = {"transcript_id": acc}

    # Data from transcript record (Phase 1 / Step 1)
    if acc in all_tracking:
        row.update({k: v for k, v in all_tracking[acc].items() if k != "transcript_id"})

    # Data from genomic record (Phase 1 / Step 2)
    genomic_acc = acc_to_genomic.get(acc)
    if genomic_acc and genomic_acc in all_tracking:
        row.update(
            {
                k: v
                for k, v in all_tracking[genomic_acc].items()
                if k != "genomic_accession"
            }
        )

    row["assembly_accession_found"] = acc_to_assembly.get(acc, "")
    row["gene_ids"] = ";".join(acc_to_geneids.get(acc, []))
    row["resolution_status"] = "resolved" if acc in resolved_set else "unresolved"
    row["unresolved_reason"] = unresolved_map.get(acc, "")
    row["resolution_strategy"] = (
        "transcript_id"
        if acc in gtf_annotations
        else "gene_id" if acc in gtf_annotations_geneid else ""
    )
    debug_rows.append(row)

debug_path = out_resolved.replace(".tsv", "_debug.tsv")
pd.DataFrame(debug_rows).to_csv(debug_path, sep="\t", index=False)
log.info(f"Written → {debug_path}")

log.info("=" * 60)
log.info(f"Input accessions          : {len(accessions)}")
log.info(f"Resolved (transcript ID)  : {n_strat_a}")
log.info(f"Resolved (Gene ID)        : {n_strat_b}")
log.info(f"Still unresolved          : {len(df_unresolved)}")
if not df_unresolved.empty:
    for reason, grp in df_unresolved.groupby("reason"):
        log.info(f"  {reason:<30}: {len(grp)}")
log.info(f"Written → {out_resolved}")
log.info(f"Debug TSV → {debug_path}")
log.info("resolve_abandoned_accessions complete.")
