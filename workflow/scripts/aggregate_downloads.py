"""
scripts/aggregate_downloads.py

Aggregate per-accession download results into summary TSVs and the
.assemblies_ready sentinel.

Reads per-accession status sentinel files (.download_done) and parses their
log files to extract detailed failure modes. Always produces all outputs so
the pipeline can proceed to the extraction and report stages even when some
downloads fail.

Snakemake interface:
    snakemake.input.status_files   — list of {CACHE}/{accession}/.download_done
    snakemake.input.log_files      — list of {LOGS}/download_assembly/{accession}.log
    snakemake.input.accession_list — results/assembly_accessions.txt
    snakemake.input.resolved       — results/ncbi_chromosome_resolved.tsv
    snakemake.output.done          — results/.assemblies_ready
    snakemake.output.downloaded    — results/downloaded_assemblies.tsv
    snakemake.output.unresolved    — results/unresolved_assemblies.tsv
    snakemake.log[0]

Failure modes reported in unresolved_assemblies.tsv:
    reason column:
        "not_resolvable_by_download_assemblies"  — non-GCF/GCA accession
        "failed: ..."                             — download failed (from status file)
    fail_detail column:
        most specific error extracted from the log file
"""

import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

log = get_logger("aggregate_downloads", snakemake.log[0])

status_files = list(snakemake.input.status_files)
log_files = list(snakemake.input.log_files)
accession_list_path = snakemake.input.accession_list
resolved_path = snakemake.input.resolved
out_done = Path(snakemake.output.done)
out_downloaded = snakemake.output.downloaded
out_unresolved = snakemake.output.unresolved


# ── Log parser ────────────────────────────────────────────────
# Patterns in priority order — first match wins.
_FAIL_PATTERNS = [
    (r"Permanent failure \(HTTP (\d+)\)[^\n]*not retrying[^\n]*", "HTTP {1} — not retrying"),
    (r"Permanent failure \(HTTP (\d+)\)",                          "HTTP {1} permanent error"),
    (r"HTTP (\d{3})",                                              "HTTP {1} error"),
    (r"Decompression failed: ([^\n]+)",                            "decompression failed: {1}"),
    (r"FAILED \(exit (\d+)\): ([^\n]+)",                           "samtools exit {1}: {2}"),
    (r"Could not determine [^\n]+URL[^\n]*",                       None),   # use full match
    (r"FTP directory listing[^\n]*contained no matching folder",   None),
    (r"All download attempts failed",                              None),
    (r"log unreadable",                                            None),
]


def parse_log_failure(log_path: str) -> str:
    """
    Return the most specific error description found in the log file.
    Falls back to the last ERROR-level line, then to "unknown error".
    """
    try:
        text = Path(log_path).read_text(errors="replace")
    except OSError as exc:
        return f"log unreadable: {exc}"

    for pattern, template in _FAIL_PATTERNS:
        m = re.search(pattern, text)
        if m:
            if template is None:
                return m.group(0).strip()
            # Fill in captured groups ({1}, {2}, ...)
            result = template
            for i, g in enumerate(m.groups(), start=1):
                result = result.replace(f"{{{i}}}", g or "")
            return result

    # Fallback: last ERROR line in the file
    for line in reversed(text.splitlines()):
        if "| ERROR    |" in line:
            # Strip the timestamp/level/name prefix, keep just the message
            parts = line.split(" | ", maxsplit=3)
            return parts[-1].strip() if len(parts) == 4 else line.strip()

    return "unknown error"


# ── Build accession → log path index ─────────────────────────
log_by_accession: dict[str, str] = {}
for lf in log_files:
    acc = Path(lf).stem   # e.g. logs/download_assembly/GCF_000001405.40.log → stem
    log_by_accession[acc] = lf

# ── Parse status files ────────────────────────────────────────
ok_accessions: set[str] = set()
failed_accessions: dict[str, str] = {}   # accession → failure detail

for sf in status_files:
    acc = Path(sf).parent.name   # resources/cache/GCF_.../  .download_done
    try:
        content = Path(sf).read_text(errors="replace").strip()
    except OSError as exc:
        content = f"failed: status file unreadable: {exc}"

    if content == "ok":
        ok_accessions.add(acc)
    else:
        # "failed: ..." from status file; augment with log details
        log_detail = parse_log_failure(log_by_accession[acc]) if acc in log_by_accession else "no log"
        failed_accessions[acc] = log_detail

log.info(f"Status files parsed       : {len(status_files)}")
log.info(f"Successfully downloaded   : {len(ok_accessions)}")
log.info(f"Download failures         : {len(failed_accessions)}")
if failed_accessions:
    for acc, detail in sorted(failed_accessions.items()):
        log.warning(f"  FAILED {acc}: {detail}")

# ── Load GCF/GCA accession list ───────────────────────────────
with open(accession_list_path) as fh:
    assembly_accessions = {line.strip() for line in fh if line.strip()}

# ── Load resolved TSV ─────────────────────────────────────────
df = pd.read_csv(resolved_path, sep="\t")

# ── Build downloaded_df ───────────────────────────────────────
downloaded_df = (
    df[df["assembly_accession"].isin(ok_accessions)]
    .drop_duplicates(subset="assembly_accession")
)

# ── Build unresolved_df ───────────────────────────────────────
# Category 1: non-GCF/GCA accessions (not in assembly_accessions.txt)
not_resolvable = df[~df["assembly_accession"].isin(assembly_accessions)].copy()
not_resolvable["reason"] = "not_resolvable_by_download_assemblies"
not_resolvable["fail_detail"] = ""
not_resolvable = not_resolvable.drop_duplicates(subset="assembly_accession")

# Category 2: GCF/GCA accessions that failed to download
failed_rows = []
for acc, detail in failed_accessions.items():
    rows = df[df["assembly_accession"] == acc]
    if rows.empty:
        # Accession was in assembly_accessions.txt but missing from resolved TSV
        row = pd.Series({"assembly_accession": acc})
        row["reason"] = "failed: download failed"
        row["fail_detail"] = detail
        failed_rows.append(row.to_frame().T)
    else:
        chunk = rows.drop_duplicates(subset="assembly_accession").copy()
        chunk["reason"] = "failed: download failed"
        chunk["fail_detail"] = detail
        failed_rows.append(chunk)

failed_df = pd.concat(failed_rows, ignore_index=True) if failed_rows else pd.DataFrame()

unresolved_df = pd.concat([not_resolvable, failed_df], ignore_index=True)

# ── Write outputs ─────────────────────────────────────────────
log.info(f"Writing {len(downloaded_df)} downloaded rows → {out_downloaded}")
downloaded_df.to_csv(out_downloaded, sep="\t", index=False)

log.info(f"Writing {len(unresolved_df)} unresolved rows → {out_unresolved}")
unresolved_df.to_csv(out_unresolved, sep="\t", index=False)

out_done.parent.mkdir(parents=True, exist_ok=True)
with open(out_done, "w") as fh:
    fh.write("assemblies_ready\n")
    fh.write(f"downloaded={len(downloaded_df)}\n")
    fh.write(f"unresolved_non_ncbi={len(not_resolvable)}\n")
    fh.write(f"download_failed={len(failed_accessions)}\n")

log.info(f"Sentinel written: {out_done}")
log.info("=" * 60)
log.info(f"Downloaded (ready for extraction) : {len(downloaded_df)}")
log.info(f"Not resolvable (non-GCF/GCA)      : {len(not_resolvable)}")
log.info(f"Download failures                 : {len(failed_accessions)}")
