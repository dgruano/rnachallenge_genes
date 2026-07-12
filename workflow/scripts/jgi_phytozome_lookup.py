"""
scripts/jgi_phytozome_lookup.py
Look up Phytozome annotation (gene GFF3) + sequence (assembly FASTA) via JGI
=============================================================================
Given a Phytozome genome_id, query the JGI file-list API and deterministically
pick the best gene-annotation GFF3 and the best assembly FASTA. Prints the
file_name, _id (needed to build a download URL), file_status, and download URL.

Usage:
  python workflow/scripts/jgi_phytozome_lookup.py --genome-id 206
  python workflow/scripts/jgi_phytozome_lookup.py --genome-id 206 --json
"""

import argparse
import json
import sys
from pathlib import Path

import requests

API_URL = "https://files.jgi.doe.gov/phytozome_file_list/"
DOWNLOAD_URL = "https://files-download.jgi.doe.gov/download_files/{_id}/"
RESTORE_URL = "https://files.jgi.doe.gov/request_archived_files/"
MAX_PAGE_SIZE = 50  # API returns HTTP 400 if x > 50 (max_page_size)


def load_token(repo_root):
    """Read JGI_SESSION_TOKEN (or PHYTOZOME_BEARER) from .env at repo root."""
    env_path = Path(repo_root) / ".env"
    if not env_path.exists():
        return None
    values = {}
    for line in env_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, val = line.partition("=")
        values[key.strip()] = val.strip().strip("'\"")
    token = values.get("JGI_SESSION_TOKEN") or values.get("PHYTOZOME_BEARER")
    if not token:
        return None
    return token if token.startswith("Bearer ") else "Bearer " + token


def fetch_files(genome_id, token, per_page=MAX_PAGE_SIZE):
    """Fetch the raw list of file dicts for a genome_id from the JGI API."""
    per_page = min(per_page, MAX_PAGE_SIZE)  # API 400s above 50
    params = [
        ("api_version", "2"),
        ("a", "false"),
        ("h", "false"),
        ("d", "asc"),
        ("p", "1"),
        ("x", str(per_page)),
        ("t", "simple"),
        ("genome_id", str(genome_id)),
    ]
    resp = requests.get(API_URL, params=params, headers={"Authorization": token})
    resp.raise_for_status()
    data = resp.json()
    organisms = data.get("organisms") or []
    if not organisms:
        return None
    return organisms[0].get("files", [])


def _file_type(f):
    return (f.get("metadata") or {}).get("type", "")


def select_files(files, prefer_name=None):
    """Pure selection: pick the best gene GFF3 and best assembly FASTA.

    Returns a dict with keys: annotation, sequence (each a file dict or None),
    plus annotation_candidates and sequence_candidates (ranked lists).

    ``prefer_name`` pins an exact portal file_name (case-insensitive) to the
    front of the ranking when present, so a manifest can name the precise file
    without breaking the default heuristic when the name is absent.
    """
    files = files or []

    annotation = []
    sequence = []
    for f in files:
        ftype = _file_type(f)
        name = (f.get("file_name") or "").lower()
        if ftype == "annotation/gene" and name.endswith(".gff3.gz"):
            annotation.append(f)
        elif ftype == "assembly" and name.endswith(".fa.gz"):
            sequence.append(f)

    def ann_key(f):
        name = (f.get("file_name") or "").lower()
        # Prefer gene_exons over plain gene; prefer versioned (contains "_v").
        return (
            "gene_exons" in name,  # True sorts after False -> reverse below
            "_v" in name,
        )

    def seq_key(f):
        name = (f.get("file_name") or "").lower()
        # Deprioritize masked assemblies.
        return "masked" not in name

    annotation.sort(key=ann_key, reverse=True)
    sequence.sort(key=seq_key, reverse=True)

    if prefer_name:
        pn = prefer_name.strip().lower()
        # Stable sort: exact-name matches (key False) float to front,
        # existing ranking preserved among the rest.
        annotation.sort(key=lambda f: (f.get("file_name") or "").lower() != pn)
        sequence.sort(key=lambda f: (f.get("file_name") or "").lower() != pn)

    return {
        "annotation": annotation[0] if annotation else None,
        "sequence": sequence[0] if sequence else None,
        "annotation_candidates": annotation,
        "sequence_candidates": sequence,
    }


def describe(f):
    """Build a small dict with the fields a caller needs to download."""
    if not f:
        return None
    return {
        "file_name": f.get("file_name"),
        "_id": f.get("_id"),
        "file_status": f.get("file_status"),
        "download_url": DOWNLOAD_URL.format(_id=f.get("_id")),
    }


def request_restore(ids, token, send_mail=False):
    """Ask JGI to restore PURGED files from tape to disk.

    Restoration can take up to 24h; JGI no-ops if the files are already on
    disk, so calling this repeatedly is safe. Returns the parsed JSON response
    (carries a request id for polling request_archived_files/requests/<id>).
    """
    if isinstance(ids, str):
        ids = [ids]
    resp = requests.post(
        RESTORE_URL,
        headers={"Authorization": token, "Content-Type": "application/json"},
        json={"ids": list(ids), "send_mail": bool(send_mail), "api_version": "2"},
    )
    resp.raise_for_status()
    return resp.json()


def resolve_annotation(genome_id, token, prefer_name=None, per_page=MAX_PAGE_SIZE):
    """Resolve a genome_id to its downloadable gene-annotation GFF3.

    Returns describe()'s dict (file_name, _id, file_status, download_url) or
    None. This is the piece the download rule needs: the numeric portal
    file_id stored in older manifests is NOT a JGI download id — only the
    Mongo _id resolved here builds a working download_files/ URL.
    """
    files = fetch_files(genome_id, token, per_page=per_page)
    if not files:
        return None
    return describe(select_files(files, prefer_name=prefer_name)["annotation"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-id", required=True, help="Phytozome genome_id")
    parser.add_argument(
        "--json", action="store_true", help="Emit machine-readable JSON"
    )
    parser.add_argument(
        "--per-page",
        type=int,
        default=MAX_PAGE_SIZE,
        help="API page size (max/default 50; higher is clamped)",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[2]
    token = load_token(repo_root)
    if not token:
        sys.exit(
            "ERROR: no JGI token found. Set JGI_SESSION_TOKEN or "
            "PHYTOZOME_BEARER in .env at the repo root."
        )

    files = fetch_files(args.genome_id, token, per_page=args.per_page)
    if files is None:
        sys.exit(f"ERROR: API returned no organisms for genome_id={args.genome_id}")

    picked = select_files(files)
    result = {
        "genome_id": args.genome_id,
        "annotation": describe(picked["annotation"]),
        "sequence": describe(picked["sequence"]),
    }

    if args.json:
        print(json.dumps(result, indent=2))
        return

    for label in ("annotation", "sequence"):
        info = result[label]
        print(f"[{label}]")
        if not info:
            print("  (none found)")
            continue
        print(f"  file_name   : {info['file_name']}")
        print(f"  _id         : {info['_id']}")
        print(f"  file_status : {info['file_status']}")
        print(f"  download_url: {info['download_url']}")
        if info["file_status"] == "PURGED":
            print("  NOTE: PURGED -> on cold storage, request restore before download")


if __name__ == "__main__":
    main()
