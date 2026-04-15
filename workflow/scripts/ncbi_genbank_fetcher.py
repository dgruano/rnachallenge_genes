"""
scripts/ncbi_genbank_fetcher.py
NCBI nuccore batch retrieval via EPost → EFetch → GenBank parsing
=================================================================
Retrieves GenBank records from NCBI's nuccore database for a list of
accession numbers, including suppressed / removed records, using
Biopython's Bio.Entrez throughout (no direct HTTP calls).

EPost → EFetch WebEnv pattern
------------------------------
1. Upload all IDs to NCBI's server-side history via epost.
2. Paginate efetch with batch_size=500 using the returned WebEnv + QueryKey.
3. Parse each batch with SeqIO.parse(handle, "genbank").
4. Log any accessions absent from the response as WITHDRAWN.

Output columns
--------------
accession | organism | description | sequence | gene_id | gene_symbol | note

Usage (standalone)
------------------
    python ncbi_genbank_fetcher.py accessions.txt out.tsv \\
        --email you@example.com [--api-key KEY] [--batch-size 500]

Usage (from another script / Snakemake helper)
----------------------------------------------
    from ncbi_genbank_fetcher import NCBIGenBankFetcher
    fetcher = NCBIGenBankFetcher(email="you@example.com", api_key="KEY")
    df = fetcher.fetch(accessions)
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
import urllib.error
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

# ── Constants ────────────────────────────────────────────────────────────────

_FETCH_DB = "nucleotide"       # efetch accepts "nucleotide" or "nuccore"
_POST_DB = "nuccore"           # epost uses "nuccore"
_GENE_FEATURE_TYPES = {"gene", "mRNA", "CDS"}
_OUTPUT_COLS = [
    "accession",
    "organism",
    "description",
    "sequence",
    "gene_id",
    "gene_symbol",
    "note",
]

# ── Helpers ───────────────────────────────────────────────────────────────────


def _get_logger(name: str, log_path: Path | str | None = None) -> logging.Logger:
    """Return a logger writing to stderr and, optionally, a file."""
    fmt = logging.Formatter(
        fmt="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    ch = logging.StreamHandler(sys.stderr)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    if log_path is not None:
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(str(log_path), mode="w", encoding="utf-8")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)

    return logger


def _extract_record_fields(record: SeqRecord) -> dict[str, Any]:
    """Extract the required fields from a single SeqRecord."""
    gene_id = ""
    gene_symbol = ""

    for feature in record.features:
        if feature.type not in _GENE_FEATURE_TYPES:
            continue
        # Gene symbol
        if not gene_symbol:
            syms = feature.qualifiers.get("gene", [])
            if syms:
                gene_symbol = syms[0]
        # GeneID from db_xref
        if not gene_id:
            for xref in feature.qualifiers.get("db_xref", []):
                if xref.startswith("GeneID:"):
                    gene_id = xref.split("GeneID:", 1)[1].strip()
                    break
        if gene_id and gene_symbol:
            break

    return {
        "accession": record.id,
        "organism": record.annotations.get("organism", ""),
        "description": record.description,
        "sequence": str(record.seq),
        "gene_id": gene_id,
        "gene_symbol": gene_symbol,
        "note": record.annotations.get("comment", ""),
    }


# ── Main class ────────────────────────────────────────────────────────────────


class NCBIGenBankFetcher:
    """
    Retrieve GenBank records from NCBI nuccore using the EPost → EFetch pattern.

    Parameters
    ----------
    email       : Registered e-mail to comply with NCBI's usage policy.
    api_key     : Optional NCBI API key (raises rate limit 3 → 10 req/s).
    batch_size  : Records per EFetch call (default 500; NCBI max is 10 000).
    log         : Pre-configured logger; if None a new one is created.
    max_retries : Per-batch HTTP retry attempts before falling back to per-record.
    retry_wait  : Base sleep (seconds) between retries; multiplied by attempt.
    """

    def __init__(
        self,
        email: str,
        api_key: str | None = None,
        batch_size: int = 500,
        log: logging.Logger | None = None,
        max_retries: int = 3,
        retry_wait: float = 5.0,
    ) -> None:
        self.batch_size = batch_size
        self.max_retries = max_retries
        self.retry_wait = retry_wait
        self.log = log or _get_logger("ncbi_genbank_fetcher")

        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

    # ── Public API ────────────────────────────────────────────────────────────

    def fetch(self, accessions: list[str]) -> pd.DataFrame:
        """
        Fetch GenBank records for *accessions* and return a DataFrame.

        Missing accessions (suppressed, withdrawn, never uploaded) are
        represented as rows with ``None`` values and ``gene_id`` == "WITHDRAWN".

        Parameters
        ----------
        accessions : List of accession strings (versioned or base).

        Returns
        -------
        pandas.DataFrame with columns: accession | organism | description |
        sequence | gene_id | gene_symbol | note
        """
        if not accessions:
            self.log.warning("fetch() called with an empty accession list")
            return pd.DataFrame(columns=_OUTPUT_COLS)

        self.log.info(f"Uploading {len(accessions)} accessions via epost …")
        webenv, query_key, matched_count = self._epost(accessions)

        rows = self._efetch_all(webenv, query_key, total=matched_count)
        returned_accessions = {row["accession"] for row in rows}

        # Detect missing / withdrawn accessions
        for acc in accessions:
            if acc not in returned_accessions:
                # Also try without version suffix
                base = acc.split(".")[0]
                if base not in returned_accessions:
                    self.log.warning(f"WITHDRAWN / not returned: {acc}")
                    rows.append(self._withdrawn_row(acc))

        df = pd.DataFrame(rows, columns=_OUTPUT_COLS)
        self.log.info(
            f"Fetch complete — {len(df)} rows "
            f"({df['gene_id'].eq('WITHDRAWN').sum()} withdrawn/missing)"
        )
        return df

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _epost(self, accessions: list[str]) -> tuple[str, str, int]:
        """Upload accession list to NCBI history; return (WebEnv, query_key, matched_count)."""
        handle = Entrez.epost(db=_POST_DB, id=",".join(accessions))
        result = Entrez.read(handle)
        handle.close()
        webenv: str = result["WebEnv"]
        query_key: str = result["QueryKey"]
        # Resolve actual count of matched records via esearch on the stored history.
        # This is the only reliable way to know how many records NCBI matched,
        # since withdrawn/suppressed IDs silently reduce the history set size.
        count_handle = Entrez.esearch(
            db=_POST_DB, term="", usehistory="y",
            webenv=webenv, query_key=query_key, retmax=0,
        )
        count_res = Entrez.read(count_handle)
        count_handle.close()
        matched: int = int(count_res.get("Count", len(accessions)))
        self.log.info(
            f"epost matched {matched}/{len(accessions)} IDs in nuccore "
            f"({len(accessions) - matched} not found / withdrawn)"
        )
        self.log.debug(f"epost → WebEnv={webenv[:30]}… QueryKey={query_key}")
        return webenv, query_key, matched

    def _efetch_all(
        self, webenv: str, query_key: str, total: int
    ) -> list[dict[str, Any]]:
        """Page through all records using the stored history set."""
        rows: list[dict[str, Any]] = []
        for start in range(0, total, self.batch_size):
            batch_rows = self._efetch_batch(
                webenv, query_key, retstart=start, retmax=self.batch_size
            )
            rows.extend(batch_rows)
            self.log.info(
                f"  Fetched {min(start + self.batch_size, total)}/{total} records"
            )
        return rows

    def _efetch_batch(
        self,
        webenv: str,
        query_key: str,
        retstart: int,
        retmax: int,
    ) -> list[dict[str, Any]]:
        """Fetch one page; fall back to per-record on HTTP error."""
        for attempt in range(1, self.max_retries + 1):
            try:
                handle = Entrez.efetch(
                    db=_FETCH_DB,
                    rettype="gb",
                    retmode="text",
                    webenv=webenv,
                    query_key=query_key,
                    retstart=retstart,
                    retmax=retmax,
                )
                records = list(SeqIO.parse(handle, "genbank"))
                handle.close()
                return [_extract_record_fields(r) for r in records]
            except urllib.error.HTTPError as exc:
                self.log.warning(
                    f"  EFetch HTTP {exc.code} on batch retstart={retstart} "
                    f"(attempt {attempt}/{self.max_retries}): {exc.reason}"
                )
                if attempt < self.max_retries:
                    time.sleep(self.retry_wait * attempt)

        # All retries exhausted — fall back to per-record fetching via direct IDs
        self.log.warning(
            f"  Batch retstart={retstart} failed after {self.max_retries} retries; "
            "falling back to per-record fetch"
        )
        return self._efetch_per_record_fallback(webenv, query_key, retstart, retmax)

    def _efetch_per_record_fallback(
        self,
        webenv: str,
        query_key: str,
        retstart: int,
        retmax: int,
    ) -> list[dict[str, Any]]:
        """Fetch records one-by-one to isolate problematic accessions."""
        rows: list[dict[str, Any]] = []
        for offset in range(retstart, retstart + retmax):
            try:
                handle = Entrez.efetch(
                    db=_FETCH_DB,
                    rettype="gb",
                    retmode="text",
                    webenv=webenv,
                    query_key=query_key,
                    retstart=offset,
                    retmax=1,
                )
                records = list(SeqIO.parse(handle, "genbank"))
                handle.close()
                if records:
                    rows.append(_extract_record_fields(records[0]))
                else:
                    self.log.warning(f"  No record returned at history offset {offset}")
            except urllib.error.HTTPError as exc:
                self.log.error(
                    f"  Per-record fetch failed at history offset {offset}: "
                    f"HTTP {exc.code} — {exc.reason}"
                )
            time.sleep(0.12)
        return rows

    @staticmethod
    def _withdrawn_row(accession: str) -> dict[str, Any]:
        return {
            "accession": accession,
            "organism": None,
            "description": None,
            "sequence": None,
            "gene_id": "WITHDRAWN",
            "gene_symbol": None,
            "note": None,
        }


# ── CLI entry-point ───────────────────────────────────────────────────────────


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Fetch GenBank records from NCBI nuccore (EPost → EFetch).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "accessions_file",
        type=Path,
        help="Plain-text file with one accession per line.",
    )
    p.add_argument(
        "output_tsv",
        type=Path,
        help="Output TSV path.",
    )
    p.add_argument("--email", required=True, help="E-mail for NCBI Entrez.")
    p.add_argument("--api-key", default=None, help="Optional NCBI API key.")
    p.add_argument("--batch-size", type=int, default=500, help="EFetch page size.")
    p.add_argument("--log", type=Path, default=None, help="Log file path.")
    p.add_argument(
        "--max-retries", type=int, default=3, help="HTTP retry attempts per batch."
    )
    p.add_argument(
        "--retry-wait", type=float, default=5.0, help="Base retry sleep in seconds."
    )
    return p


def main(argv: list[str] | None = None) -> None:
    args = _build_parser().parse_args(argv)

    log = _get_logger("ncbi_genbank_fetcher", args.log)

    accessions = [
        line.strip()
        for line in args.accessions_file.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    log.info(f"Loaded {len(accessions)} accessions from {args.accessions_file}")

    fetcher = NCBIGenBankFetcher(
        email=args.email,
        api_key=args.api_key,
        batch_size=args.batch_size,
        log=log,
        max_retries=args.max_retries,
        retry_wait=args.retry_wait,
    )
    df = fetcher.fetch(accessions)

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_tsv, sep="\t", index=False)
    log.info(f"Written {len(df)} rows → {args.output_tsv}")


# ── Sanity check (no Snakemake needed) ───────────────────────────────────────

if __name__ == "__main__":
    # If called without arguments, run the minimal self-test; otherwise run CLI.
    if len(sys.argv) == 1:
        print("Running self-test with two known accessions …", file=sys.stderr)
        import os

        email = os.environ.get("NCBI_EMAIL", "test@example.com")
        fetcher = NCBIGenBankFetcher(email=email, batch_size=500)
        test_accs = ["NM_007294.4", "NM_000546.6"]   # BRCA1, TP53 — stable records
        df = fetcher.fetch(test_accs)
        assert len(df) == len(test_accs), f"Expected {len(test_accs)} rows, got {len(df)}"
        assert "gene_id" in df.columns
        assert "gene_symbol" in df.columns
        assert df["gene_id"].ne("WITHDRAWN").any(), "All records were marked WITHDRAWN"
        print(df[["accession", "gene_symbol", "organism"]].to_string(), file=sys.stderr)
        print("Self-test PASSED", file=sys.stderr)
    else:
        main()
