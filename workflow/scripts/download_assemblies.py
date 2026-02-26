"""
scripts/download_assemblies.py
Stage 3 — Download & Cache Genome Assemblies
=============================================
Reads the resolved TSV to find all unique (organism, assembly_accession) pairs.
For each pair:
  1. Checks if already cached (skips if so)
  2. Determines the FTP source (NCBI FTP for GCF_/GCA_, Ensembl FTP otherwise)
  3. Downloads the primary genome FASTA (chromosomes / top-level)
  4. Decompresses (.gz)
  5. Indexes with samtools faidx

Cache layout:
  resources/cache/
    <assembly_accession>/
      genome.fasta
      genome.fasta.fai

Writes a sentinel file (.assemblies_ready) on success.
"""

import gzip
import hashlib
import shutil
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
out_sentinel = snakemake.output.done
cfg = snakemake.config

CACHE_DIR = Path(cfg["cache_dir"])
MAX_RETRIES = int(cfg.get("max_retries", 3))
RETRY_WAIT = int(cfg.get("retry_wait_seconds", 5))

NCBI_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
ENSEMBL_FTP_BASE = "https://ftp.ensembl.org/pub/current_fasta"

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


def decompress_gz(gz_path: Path, out_path: Path, label: str) -> bool:
    """Decompress a .gz file."""
    log.info(f"  [{label}] Decompressing {gz_path.name} ...")
    try:
        with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        gz_path.unlink()
        log.info(f"  [{label}] Decompressed → {out_path}")
        return True
    except Exception as exc:
        log.error(f"  [{label}] Decompression failed: {exc}")
        return False


def index_fasta(fasta_path: Path, label: str) -> bool:
    """Index a FASTA file with samtools faidx."""
    log.info(f"  [{label}] Indexing with samtools faidx ...")
    return run_cmd(["samtools", "faidx", str(fasta_path)], label)


# ── NCBI FTP URL resolution ───────────────────────────────────
def ncbi_assembly_url(accession: str) -> Optional[str]:
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
            return None
        report = reports[0]
        ftp_path = report.get("assembly_info", {}).get("assembly_stats", {})
        # Try to build the FTP URL from the accession pattern
        # GCF_000001405.40 → GCF/000/001/405/GCF_000001405.40_GRCh38.p14/
        acc_no_version = accession.split(".")[0]
        prefix = acc_no_version[0:3]  # GCF or GCA
        digits = acc_no_version[4:]  # 000001405
        d1, d2, d3 = digits[0:3], digits[3:6], digits[6:9]
        asm_name = report.get("assembly_info", {}).get("assembly_name", "")
        if not asm_name:
            return None
        full_name = f"{accession}_{asm_name}"
        url = (
            f"{NCBI_FTP_BASE}/{prefix}/{d1}/{d2}/{d3}/"
            f"{full_name}/{full_name}_genomic.fna.gz"
        )
        log.debug(f"  NCBI FTP URL: {url}")
        return url
    except Exception as exc:
        log.warning(f"  Could not resolve NCBI FTP URL for {accession}: {exc}")
        return None


# ── Ensembl FTP URL resolution ────────────────────────────────
def ensembl_assembly_url(organism: str, assembly: str) -> Optional[str]:
    """
    Build an Ensembl FTP URL for the toplevel genomic FASTA.
    organism: scientific name, e.g. 'homo sapiens'
    assembly: e.g. 'GRCh38'
    """
    # Ensembl uses species names like homo_sapiens
    species = organism.lower().replace(" ", "_")
    if not species:
        return None
    url = (
        f"{ENSEMBL_FTP_BASE}/{species}/dna/"
        f"{species.capitalize()}.{assembly}.dna.toplevel.fa.gz"
    )
    log.debug(f"  Ensembl FTP URL: {url}")
    return url


# ── Per-assembly download orchestrator ───────────────────────
def ensure_assembly(accession: str, organism: str, db_source: str) -> bool:
    """
    Ensure genome FASTA for the given assembly is available and indexed.
    Returns True if ready, False if download failed.
    """
    asm_dir = CACHE_DIR / accession
    fasta = asm_dir / "genome.fasta"
    fai = Path(str(fasta) + ".fai")
    label = accession

    if fai.exists() and fasta.exists():
        log.info(f"  [{label}] Already cached and indexed — skipping download")
        return True

    asm_dir.mkdir(parents=True, exist_ok=True)
    gz_path = asm_dir / "genome.fasta.gz"

    # Determine download URL based on accession prefix
    url: Optional[str] = None
    if accession.startswith("GCF_") or accession.startswith("GCA_"):
        url = ncbi_assembly_url(accession)
    elif accession in ("hg38", "hg19", "mm39", "mm10", "rn7", "dm6", "danRer11"):
        # UCSC assembly name — map to Ensembl/NCBI equivalent
        UCSC_TO_ENSEMBL = {
            "hg38": ("homo sapiens", "GRCh38"),
            "hg19": ("homo sapiens", "GRCh37"),
            "mm39": ("mus musculus", "GRCm39"),
            "mm10": ("mus musculus", "GRCm38"),
            "rn7": ("rattus norvegicus", "mRatBN7.2"),
            "dm6": ("drosophila melanogaster", "BDGP6.46"),
            "danRer11": ("danio rerio", "GRCz11"),
        }
        mapped = UCSC_TO_ENSEMBL.get(accession)
        if mapped:
            url = ensembl_assembly_url(*mapped)
    else:
        # Assume Ensembl assembly name
        url = ensembl_assembly_url(organism, accession)

    if url is None:
        log.error(f"  [{label}] Could not determine download URL")
        return False

    # Download compressed FASTA
    if not download_file(url, gz_path, label):
        return False

    # Decompress
    if not decompress_gz(gz_path, fasta, label):
        return False

    # Index
    if not index_fasta(fasta, label):
        return False

    return True


# ── Main ─────────────────────────────────────────────────────
log.info("Stage 3: Downloading and caching genome assemblies")

CACHE_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(input_tsv, sep="\t")

# Get unique (assembly_accession, organism, db_source) combos
asm_df = (
    df[["assembly_accession", "organism", "db_source"]]
    .dropna(subset=["assembly_accession"])
    .drop_duplicates(subset="assembly_accession")
)
log.info(f"Unique assemblies to ensure: {len(asm_df)}")
for _, row in asm_df.iterrows():
    log.info(f"  {row['assembly_accession']} | {row['organism']} | {row['db_source']}")

successes = 0
failures = []

for _, row in asm_df.iterrows():
    acc = str(row["assembly_accession"])
    org = str(row.get("organism", ""))
    src = str(row.get("db_source", ""))
    log.info(f"Processing assembly: {acc} ({org})")

    ok = ensure_assembly(acc, org, src)
    if ok:
        successes += 1
        log.info(f"  ✓ {acc} ready")
    else:
        failures.append(acc)
        log.error(f"  ✗ {acc} FAILED")

# Write sentinel
sentinel = Path(out_sentinel)
sentinel.parent.mkdir(parents=True, exist_ok=True)
with open(sentinel, "w") as fh:
    fh.write(f"assemblies_ready\nsuccesses={successes}\nfailures={len(failures)}\n")
    for f in failures:
        fh.write(f"FAILED: {f}\n")

# ── Summary ──────────────────────────────────────────────────
log.info("=" * 60)
log.info(f"Total assemblies needed      : {len(asm_df)}")
log.info(f"Successfully prepared        : {successes}")
log.info(f"Failed                       : {len(failures)}")
if failures:
    log.warning(f"Failed assemblies: {failures}")
    log.warning("Transcripts requiring these assemblies will be skipped in Stage 4")
log.info(f"Cache directory              : {CACHE_DIR}")
log.info("Stage 3 complete.")
