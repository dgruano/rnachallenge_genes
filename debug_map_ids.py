#!/usr/bin/env python3
"""
Debug version of map_ids_to_tools.py with verbose output.
Traces token matching for Zm00001d transcripts.
"""

import gzip
import hashlib
import io
import json
import re
import sys
import tarfile
import zipfile
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

# --- Accession patterns (copied from original) ---
_GI_RE = re.compile(
    r"gi[_|](\d+)[_|]ref[_|]((?:NM|NR|XM|XR|NP|XP|NG|NC|NT|NW|NZ)_\d+(?:\.\d+)?)[_|]?",
    re.IGNORECASE,
)

_ACCESSION_PATTERNS = [
    re.compile(r'\b((?:NM|NR|XM|XR|NP|XP|NG|NC|NT|NW|NZ)_\d+(?:\.\d+)?)\b', re.IGNORECASE),
    re.compile(r'\b(ENS[A-Z]*T\d{11}(?:\.\d+)?)\b', re.IGNORECASE),
    re.compile(r'\b(ENS[A-Z]*G\d{11}(?:\.\d+)?)\b', re.IGNORECASE),
    re.compile(r'\b(uc\d{3}[a-z]{3}\.\d+)\b', re.IGNORECASE),
    re.compile(r'\b(NON[A-Z]{3}[TG]\d+\.\d+)\b'),
    re.compile(r'\b(AT[0-9CM]G\d{5}(?:\.\d+)?)\b', re.IGNORECASE),
    re.compile(r'\b(Os\d{2}g\d{7})\b', re.IGNORECASE),
    re.compile(r'\b(LOC_Os\d{2}g\d{5})\b', re.IGNORECASE),
    re.compile(r'\b(Glyma\.\d{2}G\d{6}(?:\.\d+)?)\b', re.IGNORECASE),
    re.compile(r'\b(Zm\d{5}g\d{6})\b', re.IGNORECASE),
    re.compile(r'\b(GRMZM\w+)\b', re.IGNORECASE),
    re.compile(r'\b(Solyc\d{2}g\d{6}\.\d+\.\d+)\b', re.IGNORECASE),
    re.compile(r'\b(WBGene\d{8})\b', re.IGNORECASE),
    re.compile(r'\b(FBtr\d{7})\b', re.IGNORECASE),
    re.compile(r'\b(Y[A-P][LR]\d{3}[WC](?:_[A-Z])?)\b', re.IGNORECASE),
]

_FASTA_SUFFIXES = {".fa", ".fasta", ".faa", ".fna", ".fa.gz", ".fasta.gz", ".faa.gz", ".fna.gz"}


def normalise_id(raw_id: str) -> str:
    """Return the canonical accession from a FASTA header first token."""
    token = raw_id.lstrip(">").strip()
    m = _GI_RE.search(token)
    if m:
        return m.group(2)
    token = token.split()[0].split("|")[0]
    return token


def _all_ids_from_header(description: str) -> set[str]:
    """Extract every plausible accession from a full FASTA header description."""
    ids: set[str] = set()

    # Pass 1: primary normalised ID
    primary = normalise_id(description)
    ids.add(primary)
    ids.add(re.sub(r"\.\d+$", "", primary))

    # Pass 2: all pipe- and space-delimited tokens
    for token in re.split(r"[|\s]+", description.lstrip(">")):
        token = token.strip()
        if not token:
            continue
        m = _GI_RE.search(token)
        norm = m.group(2) if m else token
        if len(norm) >= 4:
            ids.add(norm)
            ids.add(re.sub(r"\.\d+$", "", norm))

    # Pass 3: regex scan of full description
    for pattern in _ACCESSION_PATTERNS:
        for hit in pattern.findall(description):
            ids.add(hit)
            ids.add(re.sub(r"\.\d+$", "", hit))

    return ids


def _extract_accessions_from_text(description: str) -> set[str]:
    accessions: set[str] = set()
    for pattern in _ACCESSION_PATTERNS:
        for hit in pattern.findall(description):
            accessions.add(hit)
    return accessions


def _parse_handle(fh):
    """Parse an open FASTA text handle."""
    n_sequences = 0
    raw_headers: set[str] = set()
    id_set: set[str] = set()

    for record in SeqIO.parse(fh, "fasta"):
        n_sequences += 1
        raw_headers.add(record.description)
        ids = _all_ids_from_header(record.description)
        id_set.update(ids)

    return n_sequences, raw_headers, id_set


def parse_fasta_full(fasta_path: Path):
    """Return lookup structures from a FASTA file (simplified)."""
    n_sequences = 0
    raw_headers: set[str] = set()
    id_set: set[str] = set()
    name = str(fasta_path)

    try:
        if name.endswith((".tar", ".tar.gz", ".tgz")):
            with tarfile.open(fasta_path, "r:*") as tf:
                for member in tf.getmembers():
                    if not member.isfile():
                        continue
                    ml = member.name.lower()
                    if not any(ml.endswith(s) for s in _FASTA_SUFFIXES):
                        continue
                    raw = tf.extractfile(member)
                    if raw is None:
                        continue
                    fh = io.TextIOWrapper(gzip.open(raw, "rt") if ml.endswith(".gz") else raw)
                    n, rh, ids = _parse_handle(fh)
                    n_sequences += n
                    raw_headers.update(rh)
                    id_set.update(ids)
        elif name.endswith(".zip"):
            with zipfile.ZipFile(fasta_path, "r") as zf:
                for member in zf.namelist():
                    ml = member.lower()
                    if not any(ml.endswith(s) for s in _FASTA_SUFFIXES):
                        continue
                    with zf.open(member) as raw:
                        fh = io.TextIOWrapper(gzip.open(raw, "rt") if ml.endswith(".gz") else raw)
                        n, rh, ids = _parse_handle(fh)
                        n_sequences += n
                        raw_headers.update(rh)
                        id_set.update(ids)
        elif name.endswith(".gz"):
            with gzip.open(name, "rt") as fh:
                n, rh, ids = _parse_handle(fh)
                n_sequences += n
                raw_headers.update(rh)
                id_set.update(ids)
        else:
            with open(name, "rt") as fh:
                n, rh, ids = _parse_handle(fh)
                n_sequences += n
                raw_headers.update(rh)
                id_set.update(ids)
    except Exception as exc:
        print(f"  ERROR parsing {fasta_path}: {exc}")

    return n_sequences, raw_headers, id_set


def main():
    config_path = Path("config/config.yaml")
    manifest_path = Path("results/tool_datasets_manifest.json")
    challenge_fasta = Path("resources/RNAChallenge.fa")

    print(f"\n{'='*70}")
    print("DEBUG: Map IDs to Tools - Token Matching Analysis")
    print(f"{'='*70}\n")

    # Load manifest
    with open(manifest_path) as fh:
        manifest = json.load(fh)

    tool_files: dict[str, list[Path]] = defaultdict(list)
    for entry in manifest:
        if entry.get("status") == "ok" and entry.get("file"):
            tool_files[entry["tool"]].append(Path(entry["file"]))

    # Build tool ID sets
    print("Building tool ID lookup tables...")
    tool_id_sets: dict[str, set[str]] = {}
    tool_raw_headers: dict[str, set[str]] = {}

    for tool, files in sorted(tool_files.items()):
        all_ids: set[str] = set()
        all_headers: set[str] = set()
        n_seqs = 0

        for fp in files:
            if fp.exists():
                n, rh, ids = parse_fasta_full(fp)
                n_seqs += n
                all_ids.update(ids)
                all_headers.update(rh)
            else:
                print(f"  [WARNING] {tool}: file not found: {fp}")

        tool_id_sets[tool] = all_ids
        tool_raw_headers[tool] = all_headers
        print(f"  [{tool}]: {n_seqs:,} sequences | {len(all_ids):,} unique IDs")

    # Parse RNAChallenge and find Zm00001d transcripts
    print(f"\nParsing RNAChallenge FASTA: {challenge_fasta}")
    zm_records = []

    with open(challenge_fasta) as fh:
        for record in SeqIO.parse(fh, "fasta"):
            if record.id.startswith("Zm"):
                zm_records.append(record)

    print(f"  Found {len(zm_records)} Zm00001d transcripts\n")

    if not zm_records:
        print("No Zm00001d transcripts found in RNAChallenge.fa")
        return

    # For each Zm00001d transcript, show token extraction and matching
    print(f"{'='*70}")
    print("DETAILED TOKEN MATCHING ANALYSIS")
    print(f"{'='*70}\n")

    for record in zm_records[:5]:  # Show first 5
        tid = normalise_id(record.id)
        all_ids = _all_ids_from_header(record.description)

        print(f"Transcript: {tid}")
        print(f"Full header: {record.description}")
        print(f"\nExtracted tokens (all_ids):")
        for token in sorted(all_ids):
            print(f"  - {token}")

        print(f"\nMatching by tool (multi_token pass):")
        for tool in sorted(tool_id_sets.keys()):
            id_set = tool_id_sets[tool]
            intersection = all_ids & id_set

            if intersection:
                print(f"  [{tool}] ✓ MATCH")
                print(f"    Matching tokens: {intersection}")
                # Show sample headers from this tool that contain these tokens
                sample_headers = list(tool_raw_headers[tool])[:3]
                print(f"    Sample tool headers:")
                for hdr in sample_headers:
                    print(f"      - {hdr[:80]}")
            else:
                print(f"  [{tool}] ✗ no match")

        print(f"\n{'-'*70}\n")


if __name__ == "__main__":
    main()
