"""
scripts/extract_sequences.py
Stage 4 — Extract Gene + Flanking Sequences
============================================
For each resolved transcript:
  1. Extends gene coordinates by upstream_bp / downstream_bp
     (strand-aware, clamped to chromosome bounds)
  2. Extracts the sequence from the cached indexed assembly
     using samtools faidx (via subprocess for speed on large FASTAs)
  3. Reverse-complements if gene is on the minus strand
  4. Writes:
       output.fasta  — one record per transcript
       output.bed    — one row per transcript (6-column BED + metadata)

FASTA header format:
  >{transcript_id} gene={gene_id} symbol={gene_symbol} organism={organism}
  assembly={assembly} loc={chrom}:{ext_start}-{ext_end}({strand})
  upstream={upstream_bp} downstream={downstream_bp}
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).parent))
from chrom_translation import load_chrom_translation, resolve_chrom_key
from logging_utils import get_logger

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("extract_sequences", snakemake.log[0])
input_tsv = snakemake.input.resolved
out_fasta = snakemake.output.fasta
out_bed = snakemake.output.bed
out_failed = snakemake.output.failed
UPSTREAM = int(snakemake.params.upstream)
DOWNSTREAM = int(snakemake.params.downstream)
CACHE_DIR = Path(snakemake.params.cache_dir)


# ── Helper: get chromosome length from .fai ───────────────────
def get_chrom_lengths(fai_path: Path) -> dict[str, int]:
    """Parse samtools faidx .fai file → {chrom: length}."""
    lengths = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                lengths[parts[0]] = int(parts[1])
    return lengths


# ── Helper: extract sequence via samtools faidx ───────────────
from typing import Optional


def faidx_extract_seq(
    fasta_path: Path, chrom: str, start: int, end: int
) -> Optional[str]:
    region = f"{chrom}:{start + 1}-{end}"
    cmd = ["samtools", "faidx", str(fasta_path), region]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.error(f"  samtools faidx failed for {region}: {result.stderr.strip()}")
        return None
    lines = result.stdout.strip().split("\n")
    seq = "".join(lines[1:]).upper()
    return seq if seq else None


# ── Main ─────────────────────────────────────────────────────
log.info("Stage 4: Extracting gene + flanking sequences")
log.info(f"Flanking window: -{UPSTREAM} bp upstream, +{DOWNSTREAM} bp downstream")

df = pd.read_csv(input_tsv, sep="\t")
log.info(f"Loaded {len(df)} resolved transcripts")

fasta_records: list[SeqRecord] = []
bed_rows: list[dict] = []
failed_rows: list[dict] = []

translation_cache: dict[str, dict[str, str]] = {}  # assembly -> {alias: refseq}

skipped_no_asm = 0
skipped_no_coord = 0
skipped_no_chrom = 0
skipped_no_seq = 0
extracted = 0

for idx, row in df.iterrows():
    tid = str(row["transcript_id"])
    gene_id = str(row["gene_id"])
    symbol = str(row["gene_symbol"])
    organism = str(row["organism"])
    assembly = str(row["assembly_accession"])
    db_source = str(row.get("db_source", "unknown"))
    chrom = str(row["chrom"])
    strand = str(row["strand"])

    label = f"{tid}/{assembly}"

    # Guard: skip rows with missing coordinates
    if pd.isna(row["start"]) or pd.isna(row["end"]):
        log.warning(f"  [{label}] Missing coordinates (start/end NaN) — skipping")
        skipped_no_coord += 1
        failed_rows.append({"transcript_id": tid, "assembly_accession": assembly, "chrom": chrom, "db_source": db_source, "fail_reason": "missing_coordinates"})
        continue

    start = int(row["start"])
    end = int(row["end"])

    # Swap if inverted (resolver bug)
    if start > end:
        start, end = end, start

    # Locate cached assembly
    fasta_path = CACHE_DIR / assembly / "genome.fasta"
    fai_path = Path(str(fasta_path) + ".fai")

    if not fasta_path.exists() or not fai_path.exists():
        log.warning(f"  [{label}] Assembly not cached — skipping")
        skipped_no_asm += 1
        failed_rows.append({"transcript_id": tid, "assembly_accession": assembly, "chrom": chrom, "db_source": db_source, "fail_reason": "assembly_not_cached"})
        continue

    # Get chromosome lengths for clamping
    chrom_lengths = get_chrom_lengths(fai_path)

    # Translate friendly chrom name → RefSeq accession for NCBI FASTAs
    # (memoized per assembly; {} for non-NCBI, falls through to chr-toggle)
    if assembly not in translation_cache:
        translation_cache[assembly] = load_chrom_translation(
            CACHE_DIR / assembly / "assembly_report.txt"
        )
    xlate = translation_cache[assembly]

    # Resolve chrom name to a .fai seqid (report map + 'chr' prefix toggle).
    chrom_key = resolve_chrom_key(chrom, xlate, chrom_lengths)

    if chrom_key is None:
        log.warning(f"  [{label}] Chromosome {chrom!r} not found in .fai — skipping")
        skipped_no_chrom += 1
        failed_rows.append({"transcript_id": tid, "assembly_accession": assembly, "chrom": chrom, "db_source": db_source, "fail_reason": "chrom_not_found"})
        continue
    if chrom_key != chrom:
        log.debug(f"  [{label}] chrom name remapped {chrom!r} → {chrom_key!r}")

    chrom_len = chrom_lengths[chrom_key]

    # Compute extended coordinates (strand-aware)
    if strand == "+":
        ext_start = max(0, start - UPSTREAM)
        ext_end = min(chrom_len, end + DOWNSTREAM)
    else:  # minus strand: upstream is toward higher coords
        ext_start = max(0, start - DOWNSTREAM)
        ext_end = min(chrom_len, end + UPSTREAM)

    log.debug(
        f"  [{label}] Gene: {chrom_key}:{start}-{end}({strand}) "
        f"→ Extended: {chrom_key}:{ext_start}-{ext_end}"
    )

    # Extract sequence
    raw_seq = faidx_extract_seq(fasta_path, chrom_key, ext_start, ext_end)
    if raw_seq is None:
        log.warning(f"  [{label}] Sequence extraction failed — skipping")
        skipped_no_seq += 1
        failed_rows.append({"transcript_id": tid, "assembly_accession": assembly, "chrom": chrom, "db_source": db_source, "fail_reason": "sequence_error"})
        continue

    # Reverse complement for minus strand
    bio_seq = Seq(raw_seq)
    if strand == "-":
        bio_seq = bio_seq.reverse_complement()
        log.debug(f"  [{label}] Applied reverse complement (minus strand)")

    # Build FASTA header
    header = (
        f"{tid} gene={gene_id} symbol={symbol} organism={organism} "
        f"assembly={assembly} loc={chrom_key}:{ext_start}-{ext_end}({strand}) "
        f"upstream={UPSTREAM} downstream={DOWNSTREAM}"
    )
    record = SeqRecord(bio_seq, id=tid, description=header.split(" ", 1)[1])
    fasta_records.append(record)

    # Build BED row (0-based half-open, standard BED)
    bed_rows.append(
        {
            "chrom": chrom_key,
            "chromStart": ext_start,
            "chromEnd": ext_end,
            "name": tid,
            "score": 0,
            "strand": strand,
            "gene_id": gene_id,
            "gene_symbol": symbol,
            "organism": organism,
            "assembly": assembly,
            "gene_start": start,
            "gene_end": end,
            "upstream_bp": UPSTREAM,
            "downstream_bp": DOWNSTREAM,
            "is_ambiguous": row.get("is_ambiguous", False),
        }
    )

    extracted += 1
    if extracted % 100 == 0:
        log.info(f"  Progress: {extracted}/{len(df)} extracted")

# ── Write FASTA ───────────────────────────────────────────────
log.info(f"Writing {len(fasta_records)} sequences to {out_fasta}")
SeqIO.write(fasta_records, out_fasta, "fasta")

# ── Write BED ─────────────────────────────────────────────────
df_bed = pd.DataFrame(bed_rows)
log.info(f"Writing {len(df_bed)} BED records to {out_bed}")
df_bed.to_csv(out_bed, sep="\t", index=False, header=True)

# ── Write failed transcripts ──────────────────────────────────
df_failed = pd.DataFrame(failed_rows)
log.info(f"Writing {len(df_failed)} failed records to {out_failed}")
df_failed.to_csv(out_failed, sep="\t", index=False)

# ── Summary ──────────────────────────────────────────────────
log.info("=" * 60)
log.info(f"Total transcripts input      : {len(df)}")
log.info(f"Successfully extracted       : {extracted}")
log.info(f"Skipped (assembly missing)   : {skipped_no_asm}")
log.info(f"Skipped (missing coords)     : {skipped_no_coord}")
log.info(f"Skipped (chrom not in .fai)  : {skipped_no_chrom}")
log.info(f"Skipped (sequence error)     : {skipped_no_seq}")
log.info(f"Written FASTA → {out_fasta}")
log.info(f"Written BED   → {out_bed}")
log.info("Stage 4 complete.")
