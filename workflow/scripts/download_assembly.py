"""
scripts/download_assembly.py

Per-accession NCBI assembly download (fan-out rule).

Always exits 0 — failures are recorded in the status sentinel so that
download_assemblies_done can run regardless of individual failures.

Snakemake interface:
    snakemake.wildcards.accession
    snakemake.input.manifest
    snakemake.output.status   — {CACHE}/{accession}/.download_done
    snakemake.log[0]
    snakemake.config: cache_dir, max_retries, retry_wait_seconds

Status file contents:
    "ok"                  — FASTA downloaded, decompressed, and indexed
    "failed: {reason}"    — human-readable failure reason
"""

import gzip
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional
from urllib.parse import unquote

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger
from ncbi_assembly_utils import ncbi_ftp_species_dir

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("download_assembly", snakemake.log[0])
accession = snakemake.wildcards.accession
manifest_path = Path(snakemake.input.manifest)
status_out = Path(snakemake.output.status)
cfg = snakemake.config

CACHE_DIR = Path(cfg["cache_dir"])
MAX_RETRIES = int(cfg.get("max_retries", 3))
RETRY_WAIT = int(cfg.get("retry_wait_seconds", 5))

CHUNK_SIZE = 8 * 1024 * 1024  # 8 MB download chunks


def _load_manifest_map(path: Path) -> dict[str, Optional[str]]:
    """Load cache_key -> fasta_url map from prepare_accession_list output."""
    if not path.exists():
        return {}
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as exc:
        log.warning(f"  Could not read manifest {path}: {exc}")
        return {}
    if df.empty or "cache_key" not in df.columns:
        return {}

    out: dict[str, Optional[str]] = {}
    fasta_col = "fasta_url" if "fasta_url" in df.columns else None
    for row in df.itertuples(index=False):
        key = str(getattr(row, "cache_key", "")).strip()
        if not key:
            continue
        url: Optional[str] = None
        if fasta_col:
            raw = getattr(row, "fasta_url", pd.NA)
            if pd.notna(raw):
                text = str(raw).strip()
                if text and text.lower() not in {"nan", "none"}:
                    url = text
        out[key] = url
    return out


MANIFEST_URLS = _load_manifest_map(manifest_path)


# ── Helpers ───────────────────────────────────────────────────
def run_cmd(cmd: list[str], label: str) -> bool:
    """Run a shell command; return True on success."""
    log.debug(f"  [{label}] Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.error(
            f"  [{label}] FAILED (exit {result.returncode}): {result.stderr.strip()}"
        )
        return False
    log.debug(f"  [{label}] OK")
    return True


def download_file(url: str, dest: Path, label: str) -> bool:
    """Stream-download url → dest. Returns True on success."""
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            log.info(f"  [{label}] Downloading (attempt {attempt}): {url}")
            with requests.get(url, stream=True, timeout=120) as r:
                r.raise_for_status()
                total = int(r.headers.get("content-length", 0))
                downloaded = 0
                with open(dest, "wb") as fh:
                    for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                        fh.write(chunk)
                        downloaded += len(chunk)
                if total and downloaded < total:
                    raise IOError(f"Incomplete download: {downloaded}/{total} bytes")
            log.info(f"  [{label}] Download complete → {dest}")
            return True
        except requests.HTTPError as exc:
            status = exc.response.status_code if exc.response is not None else None
            if dest.exists():
                dest.unlink()
            if status is not None and 400 <= status < 500:
                log.error(
                    f"  [{label}] Permanent failure (HTTP {status}) — not retrying: {url}"
                )
                return False
            log.warning(
                f"  [{label}] attempt {attempt} failed (HTTP {status}, transient): {exc}"
            )
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_WAIT * attempt)
        except Exception as exc:
            log.warning(f"  [{label}] attempt {attempt} failed (transient): {exc}")
            if dest.exists():
                dest.unlink()
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_WAIT * attempt)
    log.error(f"  [{label}] All download attempts failed")
    return False


def index_fasta(fasta_path: Path, label: str) -> bool:
    """Index a FASTA file with samtools faidx."""
    log.info(f"  [{label}] Indexing with samtools faidx ...")
    return run_cmd(["samtools", "faidx", str(fasta_path)], label)


def is_ncbi_assembly_accession(accession: str) -> bool:
    """Check if accession is a downloadable NCBI assembly accession (GCF_/GCA_)."""
    try:
        if pd.isna(accession):
            return False
    except (TypeError, ValueError):
        pass
    if not accession:
        return False
    acc_str = str(accession).strip()
    return acc_str.startswith(("GCF_", "GCA_"))


def ftp_assembly_folder(accession: str) -> Optional[str]:
    """Resolve an accession's NCBI FTP assembly-folder URL via directory listing.

    Reliable alternative to the datasets API (which rate-limits / 429s): lists
    the parent species dir and picks the ``{accession}_<name>/`` folder. Returns
    the folder URL (trailing slash) or None.
    """
    ftp_dir = ncbi_ftp_species_dir(accession)
    try:
        dir_resp = requests.get(ftp_dir, timeout=30)
        dir_resp.raise_for_status()
    except Exception as exc:
        log.warning(f"  FTP directory listing failed for {accession}: {exc}")
        return None
    match = re.search(rf'href="({re.escape(accession)}_[^/"]+)/"', dir_resp.text)
    if not match:
        log.warning(f"  No matching FTP folder for {accession} under {ftp_dir}")
        return None
    return f"{ftp_dir}{match.group(1)}/"


def ncbi_fasta_url(accession: str) -> Optional[str]:
    """
    Resolve an NCBI assembly accession (GCF_/GCA_) to the genomic FASTA URL.
    Prefers the datasets API summary endpoint; falls back to FTP listing.
    """
    api_url = (
        f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
        f"{accession}/dataset_report"
    )
    try:
        resp = requests.get(api_url, timeout=30)
        resp.raise_for_status()
        reports = resp.json().get("reports", [])
        asm_name = (
            reports[0].get("assembly_info", {}).get("assembly_name", "")
            if reports
            else ""
        )
        if asm_name:
            full_name = f"{accession}_{unquote(asm_name).replace(' ', '_')}"
            url = f"{ncbi_ftp_species_dir(accession)}{full_name}/{full_name}_genomic.fna.gz"
            log.debug(f"  NCBI FASTA URL: {url}")
            return url
    except Exception as exc:
        log.warning(f"  datasets API lookup failed for {accession}: {exc}")

    # Fallback: derive from the FTP directory listing (API empty or errored).
    folder = ftp_assembly_folder(accession)
    if not folder:
        log.warning(f"  Could not resolve NCBI FASTA URL for {accession}")
        return None
    name = folder.rstrip("/").rsplit("/", 1)[-1]
    return f"{folder}{name}_genomic.fna.gz"


def fetch_assembly_report(accession: str, asm_dir: Path, label: str) -> None:
    """Best-effort cache of NCBI *_assembly_report.txt for chrom-name translation.

    Resolves the report URL straight from the FTP directory listing — decoupled
    from the rate-limited datasets API, which previously 429'd on rerun and left
    already-cached genomes without a report (→ chrom_not_found). The report is
    the sibling of the genomic FASTA in the assembly folder. Non-fatal and
    idempotent: a missing report just means extract falls back to the chr-prefix
    toggle for this assembly.
    """
    report_out = asm_dir / "assembly_report.txt"
    if report_out.exists():
        return
    folder = ftp_assembly_folder(accession)
    if not folder:
        log.warning(f"  [{label}] Could not resolve FTP folder — skipping report")
        return
    name = folder.rstrip("/").rsplit("/", 1)[-1]
    report_url = f"{folder}{name}_assembly_report.txt"
    if not download_file(report_url, report_out, f"{label}:report"):
        log.warning(f"  [{label}] Assembly report unavailable (non-fatal)")


# ── Main logic (returns status string) ───────────────────────
def run_download(accession: str) -> str:
    """
    Attempt to download and index the assembly.
    Returns "ok" on success, "failed: {reason}" on any failure.
    """
    label = accession
    asm_dir = CACHE_DIR / accession
    fasta_out = asm_dir / "genome.fasta"
    fai_out = asm_dir / "genome.fasta.fai"
    fasta_gz = asm_dir / "genomic.fna.gz"
    direct_fasta_url = MANIFEST_URLS.get(accession)

    log.info(f"download_assembly: {accession}")

    # Case 1: already fully cached
    if fai_out.exists() and fasta_out.exists():
        log.info(f"  [{label}] Already cached and indexed — nothing to do")
        if is_ncbi_assembly_accession(accession):
            fetch_assembly_report(accession, asm_dir, label)
        return "ok"

    asm_dir.mkdir(parents=True, exist_ok=True)

    # Case 2: FASTA exists but .fai missing — skip straight to indexing
    if fasta_out.exists() and not fai_out.exists():
        log.info(f"  [{label}] FASTA exists but .fai missing — re-indexing")
        if not index_fasta(fasta_out, label):
            return "failed: samtools faidx failed on existing FASTA"
        if is_ncbi_assembly_accession(accession):
            fetch_assembly_report(accession, asm_dir, label)
        return "ok"

    # Case 3: full download required
    if direct_fasta_url:
        url = direct_fasta_url
    elif is_ncbi_assembly_accession(accession):
        url = ncbi_fasta_url(accession)
        if url is None:
            return "failed: could not determine NCBI FTP URL"
    else:
        return "failed: no fasta_url in manifest and not an NCBI accession"

    if not download_file(url, fasta_gz, label):
        # Extract the HTTP status from the log if available — best-effort
        return "failed: download failed (see log for details)"

    try:
        log.info(f"  [{label}] Decompressing {fasta_gz} → {fasta_out}")
        with gzip.open(fasta_gz, "rb") as f_in:
            with open(fasta_out, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        log.debug(f"  [{label}] Decompression complete")
    except Exception as exc:
        log.error(f"  [{label}] Decompression failed: {exc}")
        return f"failed: decompression error — {exc}"
    finally:
        if fasta_gz.exists():
            fasta_gz.unlink()

    if not index_fasta(fasta_out, label):
        return "failed: samtools faidx failed"

    if is_ncbi_assembly_accession(accession):
        fetch_assembly_report(accession, asm_dir, label)
    log.info(f"  [{label}] Ready: {fasta_out}")
    return "ok"


# ── Entry point — always write status, always exit 0 ─────────
status_out.parent.mkdir(parents=True, exist_ok=True)
result = run_download(accession)
status_out.write_text(result + "\n")

if result == "ok":
    log.info(f"  [{accession}] status: ok")
else:
    log.error(f"  [{accession}] status: {result}")
