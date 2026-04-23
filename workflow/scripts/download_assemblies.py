"""
scripts/download_assemblies.py

Stage 3 / Phase 4 — Download & Cache Genome Assemblies (Simplified)
====================================================================

Simplified implementation handling only NCBI assembly accessions (GCF_/GCA_).

Reads the resolved TSV to find all unique assembly_accession values.
For each GCF_/GCA_ accession:
  1. Checks if already cached (skips if present)
  2. Downloads from NCBI FTP → resources/cache/<accession>/genomic.gtf.gz
  3. Indexes with samtools faidx

Non-GCF_/GCA_ accessions are marked as unresolved with reason "not_resolvable_by_download_assemblies".

Cache layout:
  resources/cache/
    <assembly_accession>/
      genomic.gtf.gz
      genomic.gtf.gz.tbi (via samtools index)

Output files:
  results/downloaded_assemblies.tsv - assemblies successfully downloaded/cached
  results/unresolved_assemblies.tsv - non-GCF_/GCA_ accessions
"""

import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("download_assemblies", snakemake.log[0])
input_tsv = snakemake.input.resolved
output_downloaded = snakemake.output.downloaded
output_unresolved = snakemake.output.unresolved
out_sentinel = snakemake.output.done
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
        except Exception as exc:
            log.warning(f"  [{label}] attempt {attempt} failed: {exc}")
            if dest.exists():
                dest.unlink()
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_WAIT * attempt)
    log.error(f"  [{label}] All download attempts failed")
    return False


def index_gtf(gtf_path: Path, label: str) -> bool:
    """Index a GTF file with samtools index."""
    log.info(f"  [{label}] Indexing with samtools index ...")
    return run_cmd(["samtools", "index", str(gtf_path)], label)


# ── NCBI FTP URL resolution ───────────────────────────────────
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


def ncbi_gtf_url(accession: str) -> Optional[str]:
    """
    Resolve an NCBI assembly accession (GCF_/GCA_) to the genomic GTF URL.
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
            return None
        report = reports[0]
        asm_name = report.get("assembly_info", {}).get("assembly_name", "")
        if not asm_name:
            return None

        # Build FTP URL from accession pattern
        # GCF_000001405.40 → GCF/000/001/405/GCF_000001405.40_GRCh38.p14/
        acc_no_version = accession.split(".")[0]
        prefix = acc_no_version[0:3]  # GCF or GCA
        digits = acc_no_version[4:]   # 000001405
        d1, d2, d3 = digits[0:3], digits[3:6], digits[6:9]
        full_name = f"{accession}_{asm_name}"
        url = (
            f"{NCBI_FTP_BASE}/{prefix}/{d1}/{d2}/{d3}/"
            f"{full_name}/{full_name}_genomic.gtf.gz"
        )
        log.debug(f"  NCBI GTF URL: {url}")
        return url
    except Exception as exc:
        log.warning(f"  Could not resolve NCBI GTF URL for {accession}: {exc}")
        return None


# ── Per-assembly download orchestrator ───────────────────────
def ensure_assembly(accession: str) -> bool:
    """
    Ensure GTF for the given NCBI assembly is downloaded and indexed.
    Returns True if ready, False if download failed.
    """
    asm_dir = CACHE_DIR / accession
    gtf = asm_dir / "genomic.gtf.gz"
    tbi = asm_dir / "genomic.gtf.gz.tbi"
    label = accession

    if tbi.exists() and gtf.exists():
        log.info(f"  [{label}] Already cached and indexed — skipping download")
        return True

    asm_dir.mkdir(parents=True, exist_ok=True)

    # Determine download URL
    url = ncbi_gtf_url(accession)
    if url is None:
        log.error(f"  [{label}] Could not determine download URL")
        return False

    # Download GTF
    if not download_file(url, gtf, label):
        return False

    # Index
    if not index_gtf(gtf, label):
        return False

    return True


# ── Main ─────────────────────────────────────────────────────
log.info("Stage 3 / Phase 4: Download and cache NCBI genome assemblies (simplified)")

CACHE_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(input_tsv, sep="\t")

# Get unique assemblies, handling missing column gracefully
if "assembly_accession" not in df.columns:
    log.error("Input TSV missing 'assembly_accession' column")
    downloaded_df = pd.DataFrame()
    unresolved_df = df.copy()
    unresolved_df["reason"] = "missing_assembly_accession_column"
else:
    unique_asm = (
        df[["assembly_accession", "organism", "db_source"]]
        .dropna(subset=["assembly_accession"])
        .drop_duplicates(subset="assembly_accession")
    )

    downloaded = []
    unresolved = []

    log.info(f"Unique assemblies to process: {len(unique_asm)}")

    for _, row in unique_asm.iterrows():
        accession = str(row["assembly_accession"]).strip()

        if is_ncbi_assembly_accession(accession):
            log.info(f"Processing GCF_/GCA_ accession: {accession}")

            # Check if already cached
            cached = CACHE_DIR / accession / "genomic.gtf.gz"
            if cached.exists():
                log.info(f"  [{accession}] Already cached — skipping")
                downloaded.append(row)
            else:
                # Try to download
                ok = ensure_assembly(accession)
                if ok:
                    downloaded.append(row)
                    log.info(f"  ✓ {accession} ready")
                else:
                    log.error(f"  ✗ {accession} FAILED to download")
                    unresolved_row = row.copy()
                    unresolved_row["reason"] = "download_failed"
                    unresolved.append(unresolved_row)
        else:
            # Non-GCF_/GCA_ accession
            log.warning(
                f"Skipping non-GCF_/GCA_ accession: {accession} "
                "(not supported in simplified Phase 4)"
            )
            unresolved_row = row.copy()
            unresolved_row["reason"] = "not_resolvable_by_download_assemblies"
            unresolved.append(unresolved_row)

    downloaded_df = pd.DataFrame(downloaded) if downloaded else pd.DataFrame()
    unresolved_df = pd.DataFrame(unresolved) if unresolved else pd.DataFrame()

# Write output TSVs
log.info(f"Writing {len(downloaded_df)} downloaded assemblies to {output_downloaded}")
downloaded_df.to_csv(output_downloaded, sep="\t", index=False)

log.info(f"Writing {len(unresolved_df)} unresolved assemblies to {output_unresolved}")
unresolved_df.to_csv(output_unresolved, sep="\t", index=False)

# Write sentinel
sentinel = Path(out_sentinel)
sentinel.parent.mkdir(parents=True, exist_ok=True)
with open(sentinel, "w") as fh:
    fh.write("assemblies_ready\n")
    fh.write(f"downloaded={len(downloaded_df)}\n")
    fh.write(f"unresolved={len(unresolved_df)}\n")

# ── Summary ──────────────────────────────────────────────────
log.info("=" * 60)
log.info(f"Total unique assemblies processed : {len(unique_asm)}")
log.info(f"Successfully downloaded/cached    : {len(downloaded_df)}")
log.info(f"Unresolved (non-GCF_/GCA_)       : {len(unresolved_df)}")
log.info(f"Output files:")
log.info(f"  Downloaded: {output_downloaded}")
log.info(f"  Unresolved: {output_unresolved}")
log.info(f"Cache directory                  : {CACHE_DIR}")
log.info("Stage 3 / Phase 4 complete.")
