"""
scripts/resolve_abandoned_accessions.py
Third-pass resolver: abandoned (withdrawn/suppressed) NCBI accessions
======================================================================
Strategy
--------
Follows the notebook pipeline exactly:

1. ``efetch(db="nucleotide", rettype="gb")`` per transcript to get its
   GenBank record, then extract the parent genomic accession (NC_/NT_/NW_/AC_)
   from ``dbxrefs`` or the ``comment`` annotation field.
2. ``efetch(db="nucleotide", rettype="gb")`` per unique genomic accession to
   get the chromosome/scaffold record, then extract the assembly accession
   from ``dbxrefs`` (``Assembly:<acc>`` entries).
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
6. Build resolved / unresolved TSVs.

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
"""

import gzip
import re
import sys
import time
import urllib.request
from pathlib import Path
from typing import Any

from Bio import SeqIO

import pandas as pd
from Bio import Entrez

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake interface ──────────────────────────────────────────────────────
log = get_logger(
    "resolve_abandoned_accessions", snakemake.log[0]  # type: ignore[name-defined]
)
input_unresolved: str = snakemake.input.unresolved  # type: ignore[name-defined]
out_resolved: str = snakemake.output.resolved  # type: ignore[name-defined]
out_unresolved: str = snakemake.output.unresolved  # type: ignore[name-defined]
cfg: dict[str, Any] = snakemake.config  # type: ignore[name-defined]

Entrez.email = cfg["ncbi_email"]
Entrez.api_key = cfg["ncbi_api_key"]

CACHE_DIR = Path(cfg["cache_dir"])
MAX_RETRIES = int(cfg.get("max_retries", 3))
RETRY_WAIT = float(cfg.get("retry_wait_seconds", 5))
BATCH_SIZE = int(cfg.get("ncbi_batch_size", 50))

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


def fetch_genomic_accessions(accessions: list[str]) -> dict[str, str]:
    """Step 1: ``{transcript_acc: genomic_accession}``.

    One ``efetch(nucleotide, rettype=gb)`` per accession.  The GenBank record
    for a withdrawn/suppressed transcript still carries the parent genomic
    accession in its ``dbxrefs`` and ``comment`` fields.
    """
    results: dict[str, str] = {}
    total = len(accessions)
    for i, acc in enumerate(accessions, 1):
        if i % 50 == 0 or i == total:
            log.info(f"  Step 1: {i}/{total} transcripts fetched")

        def _fetch(a=acc):
            handle = Entrez.efetch(
                db="nucleotide", id=a, rettype="gb", retmode="text"
            )
            rec = SeqIO.read(handle, "genbank")
            handle.close()
            return rec

        try:
            gb = _retry(_fetch, f"efetch(nucleotide) {acc}")
        except RuntimeError as exc:
            log.warning(f"  {acc}: efetch failed — {exc}")
            time.sleep(0.11)
            continue

        genomic = _extract_genomic_accession_from_record(gb)
        if genomic:
            results[acc] = genomic
        else:
            log.debug(f"  {acc}: no genomic accession found in GenBank record")
        time.sleep(0.11)

    return results


def fetch_assembly_accessions(genomic_accs: list[str]) -> dict[str, str]:
    """Step 2: ``{genomic_accession: assembly_accession}``.

    Deduplicated: only unique genomic accessions are fetched.  The assembly
    accession is read from the ``Assembly:`` entry in ``dbxrefs`` of the
    genomic GenBank record, exactly as in the notebook.
    """
    unique = list(dict.fromkeys(genomic_accs))  # preserve order, deduplicate
    results: dict[str, str] = {}
    total = len(unique)
    for i, gacc in enumerate(unique, 1):
        if i % 20 == 0 or i == total:
            log.info(f"  Step 2: {i}/{total} genomic records fetched")

        def _fetch(g=gacc):
            handle = Entrez.efetch(
                db="nucleotide", id=g, rettype="gb", retmode="text"
            )
            rec = SeqIO.read(handle, "genbank")
            handle.close()
            return rec

        try:
            gb = _retry(_fetch, f"efetch(nucleotide/genomic) {gacc}")
        except RuntimeError as exc:
            log.warning(f"  {gacc}: efetch failed — {exc}")
            time.sleep(0.11)
            continue

        asm = _extract_assembly_from_genomic_record(gb)
        if asm:
            results[gacc] = asm
        else:
            log.warning(f"  {gacc}: no Assembly: entry in dbxrefs")
        time.sleep(0.11)

    return results





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
        time.sleep(0.11)

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
        time.sleep(0.11)

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
            "gtf_url": f"{ftp_path}/{prefix}_genomic.gtf.gz",
            "organism": organism,
        }

    return results


# ── Step 4: Download assembly GTF ────────────────────────────────────────────

def download_gtf(assembly_acc: str, gtf_url: str) -> Path | None:
    """
    Download the assembly GTF to cache; return path to ``.gtf.gz`` or ``None``.

    Caches under ``{CACHE_DIR}/{assembly_acc}/genomic.gtf.gz``.
    Skips the download if the file already exists.
    """
    asm_dir = CACHE_DIR / assembly_acc
    asm_dir.mkdir(parents=True, exist_ok=True)
    gz_path = asm_dir / "genomic.gtf.gz"

    if gz_path.exists():
        log.info(f"  GTF already cached: {gz_path}")
        return gz_path

    log.info(f"  Downloading GTF for {assembly_acc}: {gtf_url}")
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            urllib.request.urlretrieve(gtf_url, gz_path)
            log.info(f"  Downloaded → {gz_path} ({gz_path.stat().st_size / 1e6:.1f} MB)")
            return gz_path
        except Exception as exc:
            log.warning(
                f"  GTF download attempt {attempt}/{MAX_RETRIES} failed for {assembly_acc}: {exc}"
            )
            if gz_path.exists():
                gz_path.unlink()
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_WAIT * attempt)

    log.error(f"  All GTF download attempts failed for {assembly_acc}")
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

# ── Step 1: efetch each transcript → genomic accession ──────────────────────
log.info("Step 1: efetch each transcript GenBank record → extract genomic accession …")
log.info(f"  (This makes {len(accessions)} individual API calls — one per transcript)")
acc_to_genomic: dict[str, str] = fetch_genomic_accessions(accessions)
log.info(f"  Genomic accessions found: {len(acc_to_genomic)}/{len(accessions)}")

# ── Step 2: efetch each unique genomic record → assembly accession ───────────
unique_genomic = list(dict.fromkeys(acc_to_genomic.values()))
log.info(
    f"Step 2: efetch {len(unique_genomic)} unique genomic records → extract assembly accession …"
)
genomic_to_assembly: dict[str, str] = fetch_assembly_accessions(unique_genomic)
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
    path = download_gtf(asm, ftp_map[asm]["gtf_url"])
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

# ── Write outputs ─────────────────────────────────────────────────────────────
df_resolved = pd.DataFrame(resolved_rows, columns=RESOLVED_COLS)
df_unresolved = pd.DataFrame(unresolved_rows, columns=UNRESOLVED_COLS)

df_resolved.to_csv(out_resolved, sep="\t", index=False)
df_unresolved.to_csv(out_unresolved, sep="\t", index=False)

log.info("=" * 60)
log.info(f"Input accessions     : {len(accessions)}")
log.info(f"Resolved via GTF     : {len(df_resolved)}")
log.info(f"Still unresolved     : {len(df_unresolved)}")
if not df_unresolved.empty:
    for reason, grp in df_unresolved.groupby("reason"):
        log.info(f"  {reason:<30}: {len(grp)}")
log.info(f"Written → {out_resolved}")
log.info("resolve_abandoned_accessions complete.")
