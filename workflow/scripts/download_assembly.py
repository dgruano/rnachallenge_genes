"""
scripts/download_assembly.py

Per-accession NCBI assembly download (fan-out rule).

Always exits 0 — failures are recorded in the status sentinel so that
download_assemblies_done can run regardless of individual failures.

Snakemake interface:
    snakemake.wildcards.accession
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

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("download_assembly", snakemake.log[0])
accession = snakemake.wildcards.accession
status_out = Path(snakemake.output.status)
cfg = snakemake.config

CACHE_DIR = Path(cfg["cache_dir"])
MAX_RETRIES = int(cfg.get("max_retries", 3))
RETRY_WAIT = int(cfg.get("retry_wait_seconds", 5))

NCBI_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
CHUNK_SIZE = 8 * 1024 * 1024  # 8 MB download chunks


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


def ncbi_fasta_url(accession: str) -> Optional[str]:
    """
    Resolve an NCBI assembly accession (GCF_/GCA_) to the genomic FASTA URL.
    Uses the NCBI datasets API summary endpoint.
    """
    api_url = (
        f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
        f"{accession}/dataset_report"
    )
    try:
        resp = requests.get(api_url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        reports = data.get("reports", [])
        if not reports:
            # Fallback: query the NCBI FTP directory listing to find the assembly folder
            acc_no_version = accession.split(".")[0]
            prefix = acc_no_version[0:3]
            digits = acc_no_version[4:]
            d1, d2, d3 = digits[0:3], digits[3:6], digits[6:9]
            ftp_dir = f"{NCBI_FTP_BASE}/{prefix}/{d1}/{d2}/{d3}/"
            log.warning(
                f"  NCBI API returned no reports for {accession}; "
                f"falling back to FTP directory listing: {ftp_dir}"
            )
            try:
                dir_resp = requests.get(ftp_dir, timeout=30)
                dir_resp.raise_for_status()
                match = re.search(
                    rf'href="({re.escape(accession)}_[^/"]+)/"', dir_resp.text
                )
                if match:
                    folder = match.group(1)
                    url = f"{ftp_dir}{folder}/{folder}_genomic.fna.gz"
                    log.debug(f"  NCBI FTP fallback URL: {url}")
                    return url
                else:
                    log.warning(
                        f"  FTP directory listing for {accession} contained no matching folder"
                    )
                    return None
            except Exception as exc:
                log.warning(f"  FTP directory listing fallback failed for {accession}: {exc}")
                return None
        report = reports[0]
        asm_name = report.get("assembly_info", {}).get("assembly_name", "")
        if not asm_name:
            return None
        asm_name = unquote(asm_name).replace(" ", "_")

        acc_no_version = accession.split(".")[0]
        prefix = acc_no_version[0:3]
        digits = acc_no_version[4:]
        d1, d2, d3 = digits[0:3], digits[3:6], digits[6:9]
        full_name = f"{accession}_{asm_name}"
        url = (
            f"{NCBI_FTP_BASE}/{prefix}/{d1}/{d2}/{d3}/"
            f"{full_name}/{full_name}_genomic.fna.gz"
        )
        log.debug(f"  NCBI FASTA URL: {url}")
        return url
    except Exception as exc:
        log.warning(f"  Could not resolve NCBI FASTA URL for {accession}: {exc}")
        return None


def fetch_assembly_report(accession: str, asm_dir: Path, label: str) -> None:
    """Best-effort cache of NCBI *_assembly_report.txt for chrom-name translation.

    The report is the sibling of the genomic FASTA on the FTP folder; derive its
    URL by swapping the suffix. Non-fatal and idempotent: a missing report just
    means extract falls back to the chr-prefix toggle for this assembly.
    """
    report_out = asm_dir / "assembly_report.txt"
    if report_out.exists():
        return
    fasta_url = ncbi_fasta_url(accession)
    if not fasta_url:
        log.warning(f"  [{label}] No FASTA URL — skipping assembly report fetch")
        return
    report_url = fasta_url.replace("_genomic.fna.gz", "_assembly_report.txt")
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

    log.info(f"download_assembly: {accession}")

    # Case 1: already fully cached
    if fai_out.exists() and fasta_out.exists():
        log.info(f"  [{label}] Already cached and indexed — nothing to do")
        fetch_assembly_report(accession, asm_dir, label)
        return "ok"

    asm_dir.mkdir(parents=True, exist_ok=True)

    # Case 2: FASTA exists but .fai missing — skip straight to indexing
    if fasta_out.exists() and not fai_out.exists():
        log.info(f"  [{label}] FASTA exists but .fai missing — re-indexing")
        if not index_fasta(fasta_out, label):
            return "failed: samtools faidx failed on existing FASTA"
        fetch_assembly_report(accession, asm_dir, label)
        return "ok"

    # Case 3: full download required
    url = ncbi_fasta_url(accession)
    if url is None:
        return "failed: could not determine NCBI FTP URL"

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
