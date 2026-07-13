"""Helpers for deriving stable cache keys from FASTA URLs."""

from __future__ import annotations

import hashlib
import re
from urllib.parse import urlparse

_COMP_EXT_RE = re.compile(r"\.(gz|bgz|bz2|xz|zip)$", re.IGNORECASE)
_FASTA_EXT_RE = re.compile(r"\.(fa|fna|fasta|fas)$", re.IGNORECASE)
_UNSAFE_RE = re.compile(r"[^A-Za-z0-9._-]+")


def _strip_known_suffixes(filename: str) -> str:
    base = filename
    while True:
        new_base = _COMP_EXT_RE.sub("", base)
        if new_base == base:
            break
        base = new_base
    base = _FASTA_EXT_RE.sub("", base)
    return base


def build_cache_key_from_url(fasta_url: str) -> str:
    """Build a filesystem-safe cache key derived from the URL basename.

    The key embeds the normalized basename and a short URL hash to avoid
    collisions between similarly named FASTA files from different hosts.
    """
    if fasta_url is None:
        raise ValueError("fasta_url cannot be None")
    url = str(fasta_url).strip()
    if not url:
        raise ValueError("fasta_url cannot be empty")

    path = urlparse(url).path
    filename = path.rsplit("/", 1)[-1] if path else ""
    if not filename:
        filename = "assembly"

    stem = _strip_known_suffixes(filename)
    slug = _UNSAFE_RE.sub("_", stem).strip("._-")
    if not slug:
        slug = "assembly"

    digest = hashlib.sha1(url.encode("utf-8")).hexdigest()[:10]
    return f"url_{slug[:80]}_{digest}"
