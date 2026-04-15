"""
workflow/scripts/map_ids_to_tools.py
=====================================
Stage 0b — Map transcript IDs to source tools.

For every sequence in RNAChallenge.fa, determine which tool dataset(s) it
originally came from by cross-referencing the downloaded FASTAs.

Matching is attempted in multiple passes, each progressively more permissive.
The first pass that produces a hit is recorded as the match strategy.

  Pass 1 — exact_id
    Strict raw-header equality (RNAChallenge header == tool FASTA header).
    No tokenisation or separator splitting is applied in this pass.

    Pass 1b — exact_header_ws
        Header equality after whitespace-only canonicalisation
        (trim + collapse repeated whitespace). Separator characters are preserved.

    Pass 1c — gi_canonical
        GI accessions are canonicalised from both underscore and pipe styles and
        matched directly (e.g., gi|123|ref|XM_...| and gi_123_ref_XM_...).

    Pass 1d — version_aware
        For recognised accession patterns only, try exact accession first and then
        version-stripped accession (e.g., NM_001.3 -> NM_001).

  Pass 2 — multi_token
      All pipe- or space-separated tokens in the tool FASTA header are
      individually normalised and indexed.  Handles:
        • CPC2   AT1G10440|AT1G10440.1  → both gene and transcript IDs indexed
        • mRNN   ENSMUST…|ENSMUSG…|…   → gene ID also indexed
        • DeepCPP gi|12345|ref|NM_…|   → pipe-GI accession extracted

  Pass 3 — regex_scan
      A regex scans the full raw header of each tool FASTA for any known
      accession pattern (RefSeq, Ensembl transcript/gene, UCSC, NONCODE,
      plant IDs).  Catches accessions buried in description fields.

  Pass 4 — seq_hash
      MD5 of the upper-cased sequence.  Version- and ID-format-independent:
      the same sequence under a different accession in a different tool
      is still matched.

  Pass 5 — seq_hash_ut
      MD5 after upper-casing and U->T normalisation, so RNA/DNA alphabet
      differences are treated as equivalent.

Outputs
-------
results/tool_source_map.tsv  — one row per RNAChallenge transcript:
    transcript_id    : normalised ID (parse_ids.py convention)
    raw_header       : full FASTA header from RNAChallenge.fa
    tools            : comma-separated matched tool names
    n_tools          : count of matching tools
    primary_tool     : first matching tool (alphabetical)
    match_strategy   : comma-separated strategies that produced hits
                       (exact_id | exact_header_ws | gi_canonical |
                        version_aware | multi_token | regex_scan |
                        seq_hash | seq_hash_ut)

results/tool_source_stats.tsv — per-tool counts:
    tool, n_matched, n_sequences_loaded, n_possible_ids_loaded,
    n_hashes_loaded, n_hashes_ut_loaded, pct_challenge_matched,
    n_exact_id, n_exact_header_ws, n_gi_canonical, n_version_aware,
    n_multi_token, n_regex_scan, n_seq_hash, n_seq_hash_ut

results/tool_source_unmatched.tsv — IDs with no match in any tool
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

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("map_ids_to_tools", snakemake.log[0])
challenge_fasta = snakemake.input.challenge_fasta
manifest_json   = snakemake.input.manifest
out_map         = snakemake.output.tool_map
out_stats       = snakemake.output.tool_stats
out_unmatched   = snakemake.output.unmatched

# ── Accession patterns ────────────────────────────────────────
# GI-style (underscore or pipe separators)
_GI_RE = re.compile(
    r"gi[_|](\d+)[_|]ref[_|]((?:NM|NR|XM|XR|NP|XP|NG|NC|NT|NW|NZ)_\d+(?:\.\d+)?)[_|]?",
    re.IGNORECASE,
)

# All known accession patterns to scan for in full headers (Pass 3)
_ACCESSION_PATTERNS = [
    # RefSeq transcript / gene / predicted
    re.compile(r'\b((?:NM|NR|XM|XR|NP|XP|NG|NC|NT|NW|NZ)_\d+(?:\.\d+)?)\b', re.IGNORECASE),
    # Ensembl transcript IDs (any species)
    re.compile(r'\b(ENS[A-Z]*T\d{11}(?:\.\d+)?)\b', re.IGNORECASE),
    # Ensembl gene IDs (any species)
    re.compile(r'\b(ENS[A-Z]*G\d{11}(?:\.\d+)?)\b', re.IGNORECASE),
    # UCSC
    re.compile(r'\b(uc\d{3}[a-z]{3}\.\d+)\b', re.IGNORECASE),
    # NONCODE v5 transcript / gene
    re.compile(r'\b(NON[A-Z]{3}[TG]\d+\.\d+)\b'),
    # Arabidopsis TAIR
    re.compile(r'\b(AT[0-9CM]G\d{5}(?:\.\d+)?)\b', re.IGNORECASE),
    # Oryza sativa (RAP-DB)
    re.compile(r'\b(Os\d{2}g\d{7})\b', re.IGNORECASE),
    # Oryza sativa (MSU/TIGR)
    re.compile(r'\b(LOC_Os\d{2}g\d{5})\b', re.IGNORECASE),
    # Glycine max (Phytozome)
    re.compile(r'\b(Glyma\.\d{2}G\d{6}(?:\.\d+)?)\b', re.IGNORECASE),
    # Zea mays (MaizeGDB)
    re.compile(r'\b(Zm\d{5}g\d{6})\b', re.IGNORECASE),
    # Zea mays (legacy GRMZM)
    re.compile(r'\b(GRMZM\w+)\b', re.IGNORECASE),
    # Solanum lycopersicum (SGN)
    re.compile(r'\b(Solyc\d{2}g\d{6}\.\d+\.\d+)\b', re.IGNORECASE),
    # WormBase gene
    re.compile(r'\b(WBGene\d{8})\b', re.IGNORECASE),
    # FlyBase transcript
    re.compile(r'\b(FBtr\d{7})\b', re.IGNORECASE),
    # SGD yeast
    re.compile(r'\b(Y[A-P][LR]\d{3}[WC](?:_[A-Z])?)\b', re.IGNORECASE),
]

_FASTA_SUFFIXES = {
    ".fa",
    ".fasta",
    ".faa",
    ".fna",
    ".fa.gz",
    ".fasta.gz",
    ".faa.gz",
    ".fna.gz",
}


# ── ID helpers ────────────────────────────────────────────────

def normalise_id(raw_id: str) -> str:
    """
    Return the canonical accession from a FASTA header first token.
    Handles plain accessions, GI-underscore and GI-pipe formats.
    """
    token = raw_id.lstrip(">").strip()
    # Apply GI regex to the full raw token before splitting on pipe
    m = _GI_RE.search(token)
    if m:
        return m.group(2)
    token = token.split()[0].split("|")[0]
    return token


def _all_ids_from_header(description: str) -> set[str]:
    """
    Extract every plausible accession from a full FASTA header description.

    Pass 1 source: normalised first token (with GI handling).
    Pass 2 source: all pipe/whitespace tokens individually normalised.
    Pass 3 source: regex scan of the entire description string.
    """
    ids: set[str] = set()

    # Pass 1: primary normalised ID
    primary = normalise_id(description)
    ids.add(primary)
    ids.add(re.sub(r"\.\d+$", "", primary))  # version-stripped

    # Pass 2: all pipe- and space-delimited tokens
    for token in re.split(r"[|\s]+", description.lstrip(">")):
        token = token.strip()
        if not token:
            continue
        # Sub-GI check on each token
        m = _GI_RE.search(token)
        norm = m.group(2) if m else token
        if len(norm) >= 4:  # skip trivially short tokens
            ids.add(norm)
            ids.add(re.sub(r"\.\d+$", "", norm))

    # Pass 3: regex scan of full description
    for pattern in _ACCESSION_PATTERNS:
        for hit in pattern.findall(description):
            ids.add(hit)
            ids.add(re.sub(r"\.\d+$", "", hit))

    return ids


def _seq_hash(seq_str: str) -> str:
    return hashlib.md5(seq_str.upper().encode()).hexdigest()


def _seq_hash_ut(seq_str: str) -> str:
    # Treat RNA/DNA alphabet differences as equivalent by normalising U->T.
    return hashlib.md5(seq_str.upper().replace("U", "T").encode()).hexdigest()


def _canonicalize_header_ws(description: str) -> str:
    # Preserve non-whitespace separators; normalise only spacing.
    return " ".join(description.strip().split())


def _extract_accessions_from_text(description: str) -> set[str]:
    accessions: set[str] = set()
    for pattern in _ACCESSION_PATTERNS:
        for hit in pattern.findall(description):
            accessions.add(hit)
    return accessions


def _version_strip(acc: str) -> str:
    return re.sub(r"\.\d+$", "", acc)


def _extract_gi_accessions(description: str) -> set[str]:
    return {m.group(2) for m in _GI_RE.finditer(description)}


# ── FASTA parsers ─────────────────────────────────────────────


def _parse_handle(
    fh,
) -> tuple[
    int,
    set[str],
    set[str],
    set[str],
    set[str],
    set[str],
    set[str],
    dict[str, str],
    dict[str, str],
]:
    """
    Parse an open FASTA text handle.
    Returns:
        n_sequences          — number of FASTA records parsed.
      raw_headers          — full raw FASTA header descriptions.
      ws_headers           — whitespace-canonical header descriptions.
      gi_accessions        — GI-canonical accessions extracted from headers.
      accession_exact      — recognised accessions with versions as seen.
      accession_versionless— recognised accessions with versions stripped.
        id_set               — structured IDs used by multi-token fallback.
      hash_map             — exact upper-case sequence hash map.
      hash_ut_map          — upper-case + U->T hash map.
    """
    n_sequences = 0
    raw_headers: set[str] = set()
    ws_headers: set[str] = set()
    gi_accessions: set[str] = set()
    accession_exact: set[str] = set()
    accession_versionless: set[str] = set()
    id_set:        set[str] = set()
    hash_map: dict[str, str] = {}
    hash_ut_map: dict[str, str] = {}

    for record in SeqIO.parse(fh, "fasta"):
        n_sequences += 1
        # Pass 1: keep full raw header exactly as emitted by SeqIO
        raw_headers.add(record.description)
        ws_headers.add(_canonicalize_header_ws(record.description))
        gi_accessions.update(_extract_gi_accessions(record.description))

        accs = _extract_accessions_from_text(record.description)
        accession_exact.update(accs)
        accession_versionless.update({_version_strip(a) for a in accs})

        # Pass 2+: structured accession extraction
        ids = _all_ids_from_header(record.description)
        id_set.update(ids)

        # Pass 4/5: sequence hash variants
        h = _seq_hash(str(record.seq))
        h_ut = _seq_hash_ut(str(record.seq))
        hash_map[h] = normalise_id(record.id)
        hash_ut_map[h_ut] = normalise_id(record.id)

    return (
        n_sequences,
        raw_headers,
        ws_headers,
        gi_accessions,
        accession_exact,
        accession_versionless,
        id_set,
        hash_map,
        hash_ut_map,
    )


def parse_fasta_full(
    fasta_path: Path,
) -> tuple[
    int,
    set[str],
    set[str],
    set[str],
    set[str],
    set[str],
    set[str],
    dict[str, str],
    dict[str, str],
]:
    """
    Return all lookup structures from a FASTA file.
    Supports plain, gzip-compressed, TAR/TAR.GZ/TGZ, and ZIP archives.
    """
    n_sequences = 0
    raw_headers: set[str]          = set()
    ws_headers: set[str]           = set()
    gi_accessions: set[str]        = set()
    accession_exact: set[str]      = set()
    accession_versionless: set[str] = set()
    id_set:        set[str]        = set()
    hash_map:      dict[str, str]  = {}
    hash_ut_map:   dict[str, str]  = {}
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
                    fh = io.TextIOWrapper(
                        gzip.open(raw, "rt") if ml.endswith(".gz") else raw
                    )
                    n, rh, wh, gi, ae, av, ids, hm, hum = _parse_handle(fh)
                    n_sequences += n
                    raw_headers.update(rh)
                    ws_headers.update(wh)
                    gi_accessions.update(gi)
                    accession_exact.update(ae)
                    accession_versionless.update(av)
                    id_set.update(ids)
                    hash_map.update(hm)
                    hash_ut_map.update(hum)
        elif name.endswith(".zip"):
            with zipfile.ZipFile(fasta_path, "r") as zf:
                for member in zf.namelist():
                    ml = member.lower()
                    if not any(ml.endswith(s) for s in _FASTA_SUFFIXES):
                        continue
                    with zf.open(member) as raw:
                        fh = io.TextIOWrapper(
                            gzip.open(raw, "rt") if ml.endswith(".gz") else raw
                        )
                        n, rh, wh, gi, ae, av, ids, hm, hum = _parse_handle(fh)
                        n_sequences += n
                        raw_headers.update(rh)
                        ws_headers.update(wh)
                        gi_accessions.update(gi)
                        accession_exact.update(ae)
                        accession_versionless.update(av)
                        id_set.update(ids)
                        hash_map.update(hm)
                        hash_ut_map.update(hum)
        elif name.endswith(".gz"):
            with gzip.open(name, "rt") as fh:
                n, rh, wh, gi, ae, av, ids, hm, hum = _parse_handle(fh)
                n_sequences += n
                raw_headers.update(rh)
                ws_headers.update(wh)
                gi_accessions.update(gi)
                accession_exact.update(ae)
                accession_versionless.update(av)
                id_set.update(ids)
                hash_map.update(hm)
                hash_ut_map.update(hum)
        else:
            with open(name, "rt") as fh:
                n, rh, wh, gi, ae, av, ids, hm, hum = _parse_handle(fh)
                n_sequences += n
                raw_headers.update(rh)
                ws_headers.update(wh)
                gi_accessions.update(gi)
                accession_exact.update(ae)
                accession_versionless.update(av)
                id_set.update(ids)
                hash_map.update(hm)
                hash_ut_map.update(hum)
    except Exception as exc:
        log.warning(f"  Could not parse {fasta_path}: {exc}")
    return (
        n_sequences,
        raw_headers,
        ws_headers,
        gi_accessions,
        accession_exact,
        accession_versionless,
        id_set,
        hash_map,
        hash_ut_map,
    )


# ── Load manifest ────────────────────────────────────────────
log.info("Stage 0b: Mapping transcript IDs to source tools")

with open(manifest_json) as fh:
    manifest = json.load(fh)

tool_files: dict[str, list[Path]] = defaultdict(list)
for entry in manifest:
    if entry.get("status") == "ok" and entry.get("file"):
        tool_files[entry["tool"]].append(Path(entry["file"]))

# ── Build tool lookup tables ──────────────────────────────────
log.info("Building ID + sequence-hash lookup from downloaded FASTAs ...")
tool_header_sets: dict[str, set[str]] = {}
tool_header_ws_sets: dict[str, set[str]] = {}
tool_gi_sets: dict[str, set[str]] = {}
tool_accession_exact_sets: dict[str, set[str]] = {}
tool_accession_versionless_sets: dict[str, set[str]] = {}
tool_sequence_counts: dict[str, int] = {}
tool_id_sets:   dict[str, set[str]]        = {}
tool_hash_maps: dict[str, dict[str, str]]  = {}
tool_hash_ut_maps: dict[str, dict[str, str]] = {}

for tool, files in sorted(tool_files.items()):
    all_headers: set[str] = set()
    all_headers_ws: set[str] = set()
    all_gi: set[str] = set()
    all_acc_exact: set[str] = set()
    all_acc_versionless: set[str] = set()
    n_sequences = 0
    all_ids:   set[str]        = set()
    all_hashes: dict[str, str] = {}
    all_hashes_ut: dict[str, str] = {}
    for fp in files:
        if fp.exists():
            (
                n_seq,
                header_tokens,
                header_ws_tokens,
                gi_tokens,
                acc_exact_tokens,
                acc_versionless_tokens,
                ids,
                hm,
                hm_ut,
            ) = parse_fasta_full(fp)
            n_sequences += n_seq
            all_headers.update(header_tokens)
            all_headers_ws.update(header_ws_tokens)
            all_gi.update(gi_tokens)
            all_acc_exact.update(acc_exact_tokens)
            all_acc_versionless.update(acc_versionless_tokens)
            all_ids.update(ids)
            all_hashes.update(hm)
            all_hashes_ut.update(hm_ut)
        else:
            log.warning(f"  [{tool}] File not found: {fp}")
    tool_header_sets[tool] = all_headers
    tool_header_ws_sets[tool] = all_headers_ws
    tool_gi_sets[tool] = all_gi
    tool_accession_exact_sets[tool] = all_acc_exact
    tool_accession_versionless_sets[tool] = all_acc_versionless
    tool_sequence_counts[tool] = n_sequences
    tool_id_sets[tool]   = all_ids
    tool_hash_maps[tool] = all_hashes
    tool_hash_ut_maps[tool] = all_hashes_ut
    log.info(
        f"  [{tool}]: {n_sequences:,} sequences | "
        f"{len(all_ids):,} possible IDs | "
        f"{len(all_hashes):,} seq-hashes from {len(files)} file(s)"
    )

# ── Parse RNAChallenge.fa ────────────────────────────────────
log.info(f"Parsing RNAChallenge FASTA: {challenge_fasta}")
challenge_records: list[dict] = []

with open(challenge_fasta) as fh:
    for record in SeqIO.parse(fh, "fasta"):
        acc_exact = _extract_accessions_from_text(record.description)
        challenge_records.append({
            "transcript_id": normalise_id(record.id),
            "raw_header":    record.description,
            "raw_header_ws": _canonicalize_header_ws(record.description),
            "gi_ids":        _extract_gi_accessions(record.description),
            "acc_exact":     acc_exact,
            "acc_base":      {_version_strip(a) for a in acc_exact},
            "seq_hash":      _seq_hash(str(record.seq)),
            "seq_hash_ut":   _seq_hash_ut(str(record.seq)),
            # All IDs that could represent this transcript (for multi-token / regex passes)
            "all_ids":       _all_ids_from_header(record.description),
        })

log.info(f"  {len(challenge_records):,} sequences in RNAChallenge.fa")

# ── Cross-reference ───────────────────────────────────────────
log.info("Cross-referencing IDs (multi-pass) ...")

rows: list[dict] = []
tool_strategy_counts: dict[str, dict[str, int]] = {
    t: {
        "exact_id": 0,
        "exact_header_ws": 0,
        "gi_canonical": 0,
        "version_aware": 0,
        "multi_token": 0,
        "regex_scan": 0,
        "seq_hash": 0,
        "seq_hash_ut": 0,
    }
    for t in tool_id_sets
}

for rec in challenge_records:
    tid      = rec["transcript_id"]
    raw_hdr  = rec["raw_header"]
    raw_hdr_ws = rec["raw_header_ws"]
    gi_ids = rec["gi_ids"]
    acc_exact = rec["acc_exact"]
    acc_base = rec["acc_base"]
    all_ids  = rec["all_ids"]
    seq_hash = rec["seq_hash"]
    seq_hash_ut = rec["seq_hash_ut"]

    matched_tools:     list[str] = []
    matched_strategies: list[str] = []

    for tool in sorted(tool_id_sets.keys()):
        header_set = tool_header_sets[tool]
        header_ws_set = tool_header_ws_sets[tool]
        gi_set = tool_gi_sets[tool]
        acc_exact_set = tool_accession_exact_sets[tool]
        acc_base_set = tool_accession_versionless_sets[tool]
        id_set   = tool_id_sets[tool]
        hash_map = tool_hash_maps[tool]
        hash_ut_map = tool_hash_ut_maps[tool]

        strategy = None

        # Pass 1 — strict raw-header equality
        if raw_hdr in header_set:
            strategy = "exact_id"

        # Pass 1b — whitespace-canonical header equality
        if strategy is None and raw_hdr_ws in header_ws_set:
            strategy = "exact_header_ws"

        # Pass 1c — GI-canonical direct matching
        if strategy is None and gi_ids & gi_set:
            strategy = "gi_canonical"

        # Pass 1d — version-aware recognised accession matching
        if strategy is None and (acc_exact & acc_exact_set or acc_base & acc_base_set):
            strategy = "version_aware"

        # Pass 2 — multi-token: any of the IDs extracted from the
        # RNAChallenge header appear in the tool's expanded ID set
        if strategy is None and all_ids & id_set:
            strategy = "multi_token"

        # Pass 3 — regex scan: any accession found in the
        # RNAChallenge header matches any accession in the tool
        # (already covered by all_ids above, so this pass is
        #  implicitly included in Pass 2; kept for labelling clarity)

        # Pass 4 — sequence hash
        if strategy is None and seq_hash in hash_map:
            strategy = "seq_hash"

        # Pass 5 — U/T-normalised sequence hash
        if strategy is None and seq_hash_ut in hash_ut_map:
            strategy = "seq_hash_ut"

        if strategy:
            matched_tools.append(tool)
            matched_strategies.append(strategy)
            tool_strategy_counts[tool][strategy] += 1

    rows.append({
        "transcript_id":   tid,
        "raw_header":      rec["raw_header"],
        "tools":           ",".join(matched_tools),
        "n_tools":         len(matched_tools),
        "primary_tool":    matched_tools[0] if matched_tools else "",
        "match_strategy":  ",".join(dict.fromkeys(matched_strategies)),  # unique, ordered
    })

df_map = pd.DataFrame(rows, columns=[
    "transcript_id", "raw_header", "tools", "n_tools",
    "primary_tool", "match_strategy",
])

# ── Write outputs ─────────────────────────────────────────────
df_map.to_csv(out_map, sep="\t", index=False)
log.info(f"Tool source map written → {out_map}")

df_unmatched = df_map.loc[df_map["n_tools"] == 0, ["transcript_id", "raw_header"]].copy()
df_unmatched.to_csv(out_unmatched, sep="\t", index=False)
log.info(f"Unmatched IDs written  → {out_unmatched} ({len(df_unmatched):,} records)")

# ── Per-tool stats ────────────────────────────────────────────
all_tools_in_manifest = {e["tool"] for e in manifest}
stats_rows = []

for tool in sorted(all_tools_in_manifest):
    sc      = tool_strategy_counts.get(tool, {})
    n_total = sum(sc.values())
    n_sequences = tool_sequence_counts.get(tool, 0)
    n_ids = len(tool_id_sets.get(tool, set()))
    n_hashes = len(tool_hash_maps.get(tool, {}))
    n_hashes_ut = len(tool_hash_ut_maps.get(tool, {}))
    pct     = 100.0 * n_total / len(challenge_records) if challenge_records else 0.0
    stats_rows.append({
        "tool":                   tool,
        "n_matched":              n_total,
        "n_sequences_loaded":      n_sequences,
        "n_possible_ids_loaded":   n_ids,
        "n_hashes_loaded":         n_hashes,
        "n_hashes_ut_loaded":      n_hashes_ut,
        "pct_challenge_matched":  round(pct, 2),
        "n_exact_id":             sc.get("exact_id",    0),
        "n_exact_header_ws":      sc.get("exact_header_ws", 0),
        "n_gi_canonical":         sc.get("gi_canonical", 0),
        "n_version_aware":        sc.get("version_aware", 0),
        "n_multi_token":          sc.get("multi_token", 0),
        "n_regex_scan":           sc.get("regex_scan", 0),
        "n_seq_hash":             sc.get("seq_hash",    0),
        "n_seq_hash_ut":          sc.get("seq_hash_ut", 0),
    })

df_stats = pd.DataFrame(stats_rows)
df_stats.to_csv(out_stats, sep="\t", index=False)
log.info(f"Per-tool stats written → {out_stats}")

# ── Strategy summary ──────────────────────────────────────────
strat_counts = df_map["match_strategy"].value_counts()

n_matched_any   = (df_map["n_tools"] > 0).sum()
n_matched_multi = (df_map["n_tools"] > 1).sum()
n_unmatched     = (df_map["n_tools"] == 0).sum()

log.info("=" * 60)
log.info(f"RNAChallenge sequences         : {len(df_map):,}")
log.info(f"  Matched to ≥1 tool           : {n_matched_any:,} ({100*n_matched_any/len(df_map):.1f}%)")
log.info(f"  Matched to >1 tool           : {n_matched_multi:,}")
log.info(f"  Unmatched (no tool found)    : {n_unmatched:,}")
log.info("")
log.info("Matches by primary strategy:")
for strat in (
    "exact_id",
    "exact_header_ws",
    "gi_canonical",
    "version_aware",
    "multi_token",
    "seq_hash",
    "seq_hash_ut",
):
    n = (df_map["match_strategy"].str.contains(strat, na=False)).sum()
    log.info(f"  {strat:<15} : {n:,}")
log.info("")
log.info(f"{'Tool':<15} {'Seqs':>10} {'IDs':>12} {'Hashes':>10} {'Matched':>9} {'%Chall':>7}  "
         f"{'exact':>6} {'ex_ws':>6} {'gi':>6} {'ver':>6} {'multi':>6} {'hash':>6} {'h_ut':>6}")
log.info("-" * 94)
for row in stats_rows:
    log.info(
        f"  {row['tool']:<13} {row['n_sequences_loaded']:>10,} {row['n_possible_ids_loaded']:>12,}"
        f" {row['n_hashes_loaded']:>10,} {row['n_matched']:>9,}"
        f"  {row['pct_challenge_matched']:>5.1f}%"
        f"  {row['n_exact_id']:>6,} {row['n_exact_header_ws']:>6,}"
        f" {row['n_gi_canonical']:>6,} {row['n_version_aware']:>6,}"
        f" {row['n_multi_token']:>6,} {row['n_seq_hash']:>6,} {row['n_seq_hash_ut']:>6,}"
    )
log.info("Stage 0b complete.")
