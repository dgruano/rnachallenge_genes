#!/usr/bin/env python3
"""
Test script to evaluate faster alternatives for NCBI transcript ID resolution.

Systematically tests:
1. Batch efetch requests
2. NCBI bulk catalog files (gene2refseq, gene2accession)
3. Different rettype/retmode formats
4. epost + efetch strategy
5. Performance comparisons

Usage:
    python test_ncbi_speedup.py
"""

import gzip
import os
import sqlite3
import time
import urllib.request
from collections import defaultdict
from io import BytesIO, StringIO
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from Bio import Entrez

# ── Configuration ──────────────────────────────────────────────
NCBI_EMAIL = "daniel.garciaruano@ibgc.cnrs.fr"
NCBI_API_KEY = "fdc3c3e0bf43bfebc561e789a0884b879308"

Entrez.email = NCBI_EMAIL
if NCBI_API_KEY:
    Entrez.api_key = NCBI_API_KEY

# Test dataset: sample IDs from the failed resolution logs
# Mix of suppressed records and normal records
TEST_IDS = [
    # Suppressed records (from investigation report)
    "XM_020345473.1",
    "XM_020345467.1",
    "XM_020315385.1",
    "XM_020294737.1",
    "XM_020313512.1",
    "XM_020315948.1",
    "XM_020316852.1",
    "XM_020323553.1",
    "XM_020324435.1",
    "XM_020327859.1",
    # Normal records (should resolve via standard methods)
    "XM_006826593.2",
    "XM_006836868.3",
    "XM_006838150.3",
    "XM_006840316.2",
    "XM_006844822.3",
    "XM_006846540.3",
    "XM_006852421.3",
    "XM_006852697.3",
    "XM_006855731.3",
    "XM_006858361.3",
]

# URLs for bulk NCBI files
GENE2REFSEQ_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz"
GENE2ACCESSION_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2accession.gz"

# Working directory for test files
TEST_DIR = Path(__file__).parent / "ncbi_speedup_tests"
TEST_DIR.mkdir(exist_ok=True)

# Results storage
RESULTS = {
    "batch_efetch": {},
    "bulk_files": {},
    "formats": {},
    "epost_efetch": {},
}


# ════════════════════════════════════════════════════════════════════════════
# TEST 1: Batch efetch Requests
# ════════════════════════════════════════════════════════════════════════════


def test_batch_efetch(ids: List[str], batch_size: int = 10) -> Dict:
    """
    Test if efetch can batch multiple IDs in a single request.

    Returns:
        dict with keys: success, time_elapsed, results, errors
    """
    print(f"\n{'='*80}")
    print(f"TEST 1: Batch efetch (batch_size={batch_size})")
    print(f"{'='*80}")

    results = []
    errors = []
    start_time = time.time()

    try:
        # Try batching IDs with comma separation
        id_string = ",".join(ids[:batch_size])
        print(f"Fetching batch: {id_string}")

        handle = Entrez.efetch(
            db="nucleotide", id=id_string, rettype="gb", retmode="text"
        )

        batch_text = handle.read()
        handle.close()

        # Count how many records we got back (look for LOCUS lines)
        record_count = batch_text.count("\nLOCUS ")
        print(f"✓ Successfully fetched {record_count} records")
        print(f"  Total size: {len(batch_text):,} bytes")

        # Try to parse gene information from each record
        # Split on LOCUS to separate records
        records = []
        current_record = []

        for line in batch_text.split("\n"):
            if line.startswith("LOCUS ") and current_record:
                records.append("\n".join(current_record))
                current_record = [line]
            else:
                current_record.append(line)

        if current_record:
            records.append("\n".join(current_record))

        print(f"  Split into {len(records)} separate records")

        # Parse each record for gene information
        for i, record_text in enumerate(records):
            gene_info = parse_genbank_text(record_text)
            if gene_info:
                results.append(gene_info)
                print(
                    f"  [{i+1}] {gene_info['accession']} → Gene {gene_info['gene_id']} ({gene_info['gene_symbol']})"
                )
            else:
                print(f"  [{i+1}] Failed to parse gene info")

        elapsed = time.time() - start_time

        return {
            "success": True,
            "time_elapsed": elapsed,
            "records_fetched": record_count,
            "records_parsed": len(results),
            "results": results,
            "errors": errors,
            "bytes_downloaded": len(batch_text),
            "time_per_record": elapsed / max(record_count, 1),
        }

    except Exception as e:
        elapsed = time.time() - start_time
        errors.append(str(e))
        print(f"✗ Batch efetch failed: {e}")

        return {
            "success": False,
            "time_elapsed": elapsed,
            "results": results,
            "errors": errors,
        }


def parse_genbank_text(record_text: str) -> Optional[Dict]:
    """Parse GenBank flat file text to extract gene information."""
    gene_info = {
        "accession": "",
        "gene_id": "",
        "gene_symbol": "",
        "organism": "",
        "chrom": "",
    }

    # Extract accession from LOCUS or VERSION line
    for line in record_text.split("\n"):
        if line.startswith("VERSION"):
            parts = line.split()
            if len(parts) >= 2:
                gene_info["accession"] = parts[1]
        elif line.startswith("  ORGANISM"):
            gene_info["organism"] = line.split("ORGANISM")[1].strip()

    # Parse FEATURES for gene information
    in_gene_feature = False
    in_source_feature = False

    for line in record_text.split("\n"):
        line_stripped = line.rstrip()

        if line_stripped.startswith("     gene "):
            in_gene_feature = True
            in_source_feature = False
        elif line_stripped.startswith("     source "):
            in_source_feature = True
            in_gene_feature = False
        elif line_stripped.startswith("     ") and not line_stripped.startswith(
            "                     "
        ):
            in_gene_feature = False
            in_source_feature = False

        if in_gene_feature:
            if '/db_xref="GeneID:' in line:
                gene_info["gene_id"] = line.split("GeneID:")[1].split('"')[0]
            elif '/gene="' in line:
                gene_info["gene_symbol"] = line.split('/gene="')[1].split('"')[0]

        if in_source_feature:
            if '/chromosome="' in line:
                gene_info["chrom"] = line.split('/chromosome="')[1].split('"')[0]

    # Return None if critical fields missing
    if not gene_info["gene_id"] or not gene_info["accession"]:
        return None

    return gene_info


# ════════════════════════════════════════════════════════════════════════════
# TEST 2: NCBI Bulk Catalog Files
# ════════════════════════════════════════════════════════════════════════════


def download_if_needed(url: str, local_path: Path) -> bool:
    """Download file if it doesn't exist locally."""
    if local_path.exists():
        print(f"  Using cached file: {local_path}")
        return True

    print(f"  Downloading: {url}")
    try:
        urllib.request.urlretrieve(url, local_path)
        print(f"  ✓ Downloaded to: {local_path} ({local_path.stat().st_size:,} bytes)")
        return True
    except Exception as e:
        print(f"  ✗ Download failed: {e}")
        return False


def test_gene2refseq(ids: List[str]) -> Dict:
    """
    Test resolution using NCBI's gene2refseq.gz bulk catalog.

    File format (tab-separated):
    #tax_id  GeneID  status  RNA_nucleotide_accession.version  ...
    """
    print(f"\n{'='*80}")
    print(f"TEST 2a: gene2refseq.gz Bulk Catalog")
    print(f"{'='*80}")

    local_file = TEST_DIR / "gene2refseq.gz"

    if not download_if_needed(GENE2REFSEQ_URL, local_file):
        return {"success": False, "errors": ["Download failed"]}

    start_time = time.time()
    results = []

    # Build lookup index: accession → gene info
    print("  Building index from gene2refseq...")
    index = {}
    lines_processed = 0

    try:
        with gzip.open(local_file, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue

                lines_processed += 1
                parts = line.strip().split("\t")

                if len(parts) < 4:
                    continue

                tax_id, gene_id, status, rna_acc = parts[:4]

                # Index both with and without version
                if rna_acc and rna_acc != "-":
                    base_acc = rna_acc.split(".")[0] if "." in rna_acc else rna_acc

                    gene_data = {
                        "gene_id": gene_id,
                        "tax_id": tax_id,
                        "status": status,
                        "rna_accession": rna_acc,
                    }

                    index[rna_acc] = gene_data
                    index[base_acc] = gene_data

                if lines_processed % 1000000 == 0:
                    print(
                        f"    Processed {lines_processed:,} lines, index size: {len(index):,}"
                    )

        print(
            f"  ✓ Index built: {len(index):,} accessions from {lines_processed:,} lines"
        )

        # Lookup our test IDs
        print("\n  Looking up test IDs:")
        found = 0
        for tid in ids:
            if tid in index:
                found += 1
                gene_data = index[tid]
                results.append(
                    {
                        "transcript_id": tid,
                        "gene_id": gene_data["gene_id"],
                        "source": "gene2refseq",
                        "status": gene_data["status"],
                    }
                )
                print(
                    f"    ✓ {tid} → Gene {gene_data['gene_id']} (status: {gene_data['status']})"
                )
            else:
                # Try base ID
                base_tid = tid.split(".")[0]
                if base_tid in index and base_tid != tid:
                    found += 1
                    gene_data = index[base_tid]
                    results.append(
                        {
                            "transcript_id": tid,
                            "gene_id": gene_data["gene_id"],
                            "source": "gene2refseq",
                            "status": gene_data["status"],
                            "note": "matched_via_base_id",
                        }
                    )
                    print(f"    ✓ {tid} → Gene {gene_data['gene_id']} (via base ID)")
                else:
                    print(f"    ✗ {tid} NOT FOUND")

        elapsed = time.time() - start_time

        return {
            "success": True,
            "time_elapsed": elapsed,
            "index_size": len(index),
            "found": found,
            "total": len(ids),
            "coverage": found / len(ids) * 100,
            "results": results,
        }

    except Exception as e:
        elapsed = time.time() - start_time
        print(f"✗ Error processing gene2refseq: {e}")
        return {
            "success": False,
            "time_elapsed": elapsed,
            "errors": [str(e)],
        }


def test_gene2accession(ids: List[str]) -> Dict:
    """
    Test resolution using NCBI's gene2accession.gz bulk catalog.

    File format (tab-separated):
    #tax_id  GeneID  status  RNA_nucleotide_accession.version  ...
    Similar to gene2refseq but includes more metadata.
    """
    print(f"\n{'='*80}")
    print(f"TEST 2b: gene2accession.gz Bulk Catalog")
    print(f"{'='*80}")

    local_file = TEST_DIR / "gene2accession.gz"

    if not download_if_needed(GENE2ACCESSION_URL, local_file):
        return {"success": False, "errors": ["Download failed"]}

    start_time = time.time()
    results = []

    print("  Building index from gene2accession...")
    index = {}
    lines_processed = 0

    try:
        with gzip.open(local_file, "rt") as f:
            # File format:
            # tax_id GeneID status RNA_nucleotide_accession.version RNA_nucleotide_gi ...
            for line in f:
                if line.startswith("#"):
                    continue

                lines_processed += 1
                parts = line.strip().split("\t")

                if len(parts) < 4:
                    continue

                tax_id, gene_id, status, rna_acc = parts[:4]

                if rna_acc and rna_acc != "-":
                    base_acc = rna_acc.split(".")[0] if "." in rna_acc else rna_acc

                    gene_data = {
                        "gene_id": gene_id,
                        "tax_id": tax_id,
                        "status": status,
                        "rna_accession": rna_acc,
                    }

                    index[rna_acc] = gene_data
                    index[base_acc] = gene_data

                if lines_processed % 1000000 == 0:
                    print(
                        f"    Processed {lines_processed:,} lines, index size: {len(index):,}"
                    )

        print(
            f"  ✓ Index built: {len(index):,} accessions from {lines_processed:,} lines"
        )

        # Lookup test IDs
        print("\n  Looking up test IDs:")
        found = 0
        for tid in ids:
            if tid in index:
                found += 1
                gene_data = index[tid]
                results.append(
                    {
                        "transcript_id": tid,
                        "gene_id": gene_data["gene_id"],
                        "source": "gene2accession",
                        "status": gene_data["status"],
                    }
                )
                print(
                    f"    ✓ {tid} → Gene {gene_data['gene_id']} (status: {gene_data['status']})"
                )
            else:
                base_tid = tid.split(".")[0]
                if base_tid in index and base_tid != tid:
                    found += 1
                    gene_data = index[base_tid]
                    results.append(
                        {
                            "transcript_id": tid,
                            "gene_id": gene_data["gene_id"],
                            "source": "gene2accession",
                            "status": gene_data["status"],
                            "note": "matched_via_base_id",
                        }
                    )
                    print(f"    ✓ {tid} → Gene {gene_data['gene_id']} (via base ID)")
                else:
                    print(f"    ✗ {tid} NOT FOUND")

        elapsed = time.time() - start_time

        return {
            "success": True,
            "time_elapsed": elapsed,
            "index_size": len(index),
            "found": found,
            "total": len(ids),
            "coverage": found / len(ids) * 100,
            "results": results,
        }

    except Exception as e:
        elapsed = time.time() - start_time
        print(f"✗ Error processing gene2accession: {e}")
        return {
            "success": False,
            "time_elapsed": elapsed,
            "errors": [str(e)],
        }


# ════════════════════════════════════════════════════════════════════════════
# TEST 3: Different Entrez rettype/retmode Formats
# ════════════════════════════════════════════════════════════════════════════


def test_efetch_formats(test_id: str) -> Dict:
    """
    Compare different efetch formats for size and parsing complexity.

    Tests:
    - rettype=gb, retmode=text (GenBank flat file - current)
    - rettype=gbc, retmode=xml (GenBank XML)
    - rettype=ft (feature table)
    - rettype=native (ASN.1)
    """
    print(f"\n{'='*80}")
    print(f"TEST 3: Different efetch Formats (using {test_id})")
    print(f"{'='*80}")

    formats_to_test = [
        ("gb", "text", "GenBank flat file"),
        ("gbc", "xml", "GenBank XML"),
        ("ft", "text", "Feature table"),
        ("gb", "xml", "GenBank XML (gb)"),
    ]

    results = {}

    for rettype, retmode, description in formats_to_test:
        print(f"\n  Testing: {description} (rettype={rettype}, retmode={retmode})")

        try:
            start_time = time.time()

            handle = Entrez.efetch(
                db="nucleotide", id=test_id, rettype=rettype, retmode=retmode
            )

            data = handle.read()
            handle.close()

            elapsed = time.time() - start_time

            # Try to detect if we got gene information
            has_gene_id = False
            gene_id = ""

            if isinstance(data, bytes):
                data_str = data.decode("utf-8", errors="ignore")
            else:
                data_str = str(data)

            if "GeneID:" in data_str or "db_xref" in data_str:
                has_gene_id = True
                # Try to extract gene ID
                if "GeneID:" in data_str:
                    parts = data_str.split("GeneID:")
                    if len(parts) > 1:
                        gene_id = parts[1].split('"')[0].split()[0].split("<")[0]

            results[f"{rettype}_{retmode}"] = {
                "description": description,
                "size_bytes": len(data),
                "time_elapsed": elapsed,
                "has_gene_id": has_gene_id,
                "gene_id": gene_id,
            }

            print(f"    ✓ Size: {len(data):,} bytes")
            print(f"    ✓ Time: {elapsed:.3f}s")
            print(
                f"    ✓ Gene ID found: {has_gene_id} ({gene_id if gene_id else 'N/A'})"
            )

        except Exception as e:
            print(f"    ✗ Failed: {e}")
            results[f"{rettype}_{retmode}"] = {
                "description": description,
                "error": str(e),
            }

    return results


# ════════════════════════════════════════════════════════════════════════════
# TEST 4: epost + efetch Strategy
# ════════════════════════════════════════════════════════════════════════════


def test_epost_efetch(ids: List[str]) -> Dict:
    """
    Test using epost to upload ID list to history server,
    then bulk fetch with efetch.
    """
    print(f"\n{'='*80}")
    print(f"TEST 4: epost + efetch Strategy ({len(ids)} IDs)")
    print(f"{'='*80}")

    start_time = time.time()

    try:
        # Step 1: Upload IDs to history server
        print("  Step 1: Uploading IDs to NCBI history server...")
        id_string = ",".join(ids)

        upload_handle = Entrez.epost(db="nucleotide", id=id_string)
        upload_result = Entrez.read(upload_handle)
        upload_handle.close()

        webenv = upload_result["WebEnv"]
        query_key = upload_result["QueryKey"]

        print(f"    ✓ Uploaded {len(ids)} IDs")
        print(f"    WebEnv: {webenv[:50]}...")
        print(f"    QueryKey: {query_key}")

        # Step 2: Fetch all records using history
        print("\n  Step 2: Fetching records using history server...")

        fetch_handle = Entrez.efetch(
            db="nucleotide",
            webenv=webenv,
            query_key=query_key,
            rettype="gb",
            retmode="text",
        )

        batch_text = fetch_handle.read()
        fetch_handle.close()

        # Count records
        record_count = batch_text.count("\nLOCUS ")

        elapsed = time.time() - start_time

        print(f"    ✓ Fetched {record_count} records")
        print(f"    ✓ Total size: {len(batch_text):,} bytes")
        print(f"    ✓ Time: {elapsed:.3f}s ({elapsed/len(ids):.3f}s per ID)")

        return {
            "success": True,
            "time_elapsed": elapsed,
            "records_fetched": record_count,
            "bytes_downloaded": len(batch_text),
            "time_per_record": elapsed / max(record_count, 1),
        }

    except Exception as e:
        elapsed = time.time() - start_time
        print(f"  ✗ epost+efetch failed: {e}")

        return {
            "success": False,
            "time_elapsed": elapsed,
            "error": str(e),
        }


# ════════════════════════════════════════════════════════════════════════════
# TEST 5: Performance Comparison - Simulate 9,000 IDs
# ════════════════════════════════════════════════════════════════════════════


def estimate_performance_9k_ids():
    """
    Estimate time to resolve 9,000 IDs using each approach.
    """
    print(f"\n{'='*80}")
    print(f"PERFORMANCE ESTIMATES FOR 9,000 IDs")
    print(f"{'='*80}\n")

    TARGET_COUNT = 9000

    # Get timing data from tests
    batch_efetch = RESULTS.get("batch_efetch", {})
    bulk_files = RESULTS.get("bulk_files", {})
    epost_efetch = RESULTS.get("epost_efetch", {})

    estimates = []

    # Current approach (individual efetch)
    estimates.append(
        {
            "approach": "Current (individual efetch)",
            "time_per_id": 1.0,  # ~1 second per ID from user requirements
            "total_time_seconds": TARGET_COUNT * 1.0,
            "total_time_minutes": TARGET_COUNT * 1.0 / 60,
            "speedup": "1x (baseline)",
        }
    )

    # Batch efetch (if successful)
    if batch_efetch.get("success"):
        time_per_record = batch_efetch.get("time_per_record", 1.0)
        total_time = TARGET_COUNT * time_per_record
        speedup = (TARGET_COUNT * 1.0) / total_time
        estimates.append(
            {
                "approach": "Batch efetch (10 IDs/request)",
                "time_per_id": time_per_record,
                "total_time_seconds": total_time,
                "total_time_minutes": total_time / 60,
                "speedup": f"{speedup:.1f}x",
            }
        )

    # epost + efetch (if successful)
    if epost_efetch.get("success"):
        time_per_record = epost_efetch.get("time_per_record", 1.0)
        # epost can handle larger batches, estimate 100 IDs per batch
        num_batches = TARGET_COUNT / 100
        total_time = num_batches * (time_per_record * 100)
        speedup = (TARGET_COUNT * 1.0) / total_time
        estimates.append(
            {
                "approach": "epost + efetch (100 IDs/batch)",
                "time_per_id": time_per_record,
                "total_time_seconds": total_time,
                "total_time_minutes": total_time / 60,
                "speedup": f"{speedup:.1f}x",
            }
        )

    # Bulk files (gene2refseq or gene2accession)
    gene2refseq_result = bulk_files.get("gene2refseq", {})
    if gene2refseq_result.get("success"):
        # Download time (one-time) + instant lookups
        download_time = gene2refseq_result.get("time_elapsed", 0)
        lookup_time_per_id = 0.0001  # Virtually instant for dict lookup
        total_time = download_time + (TARGET_COUNT * lookup_time_per_id)
        speedup = (TARGET_COUNT * 1.0) / total_time
        coverage = gene2refseq_result.get("coverage", 0)

        estimates.append(
            {
                "approach": f"gene2refseq.gz (one-time download, {coverage:.1f}% coverage)",
                "time_per_id": lookup_time_per_id,
                "total_time_seconds": total_time,
                "total_time_minutes": total_time / 60,
                "speedup": f"{speedup:.0f}x",
                "note": f"Download time: {download_time:.1f}s, Coverage: {coverage:.1f}%",
            }
        )

    # Print comparison table
    print(f"{'Approach':<50} | {'Time/ID':>10} | {'Total Time':>12} | {'Speedup':>10}")
    print(f"{'-'*50}-+-{'-'*10}-+-{'-'*12}-+-{'-'*10}")

    for est in estimates:
        approach = est["approach"]
        time_per = f"{est['time_per_id']:.4f}s"
        total = (
            f"{est['total_time_minutes']:.1f} min"
            if est["total_time_minutes"] < 120
            else f"{est['total_time_minutes']/60:.1f} hr"
        )
        speedup = est["speedup"]

        print(f"{approach:<50} | {time_per:>10} | {total:>12} | {speedup:>10}")

        if "note" in est:
            print(f"  Note: {est['note']}")

    print()


# ════════════════════════════════════════════════════════════════════════════
# Main Test Runner
# ════════════════════════════════════════════════════════════════════════════


def main():
    print(f"\n{'#'*80}")
    print(f"# NCBI Transcript ID Resolution - Performance Testing")
    print(f"# Testing faster alternatives for ~9,000 transcript IDs")
    print(f"{'#'*80}\n")

    print(f"Test IDs: {len(TEST_IDS)}")
    print(f"Working directory: {TEST_DIR}")
    print(f"NCBI Email: {NCBI_EMAIL}")
    print(f"API Key configured: {'Yes' if NCBI_API_KEY else 'No'}")

    # Run tests
    try:
        # Test 1: Batch efetch
        batch_result = test_batch_efetch(TEST_IDS, batch_size=10)
        RESULTS["batch_efetch"] = batch_result
        time.sleep(1)  # Be nice to NCBI

        # Test 2a: gene2refseq
        gene2refseq_result = test_gene2refseq(TEST_IDS)
        if not RESULTS.get("bulk_files"):
            RESULTS["bulk_files"] = {}
        RESULTS["bulk_files"]["gene2refseq"] = gene2refseq_result

        # Test 2b: gene2accession
        gene2accession_result = test_gene2accession(TEST_IDS)
        RESULTS["bulk_files"]["gene2accession"] = gene2accession_result

        # Test 3: Different formats (use first test ID)
        time.sleep(1)
        formats_result = test_efetch_formats(TEST_IDS[0])
        RESULTS["formats"] = formats_result

        # Test 4: epost + efetch
        time.sleep(1)
        epost_result = test_epost_efetch(TEST_IDS)
        RESULTS["epost_efetch"] = epost_result

        # Test 5: Performance estimates
        estimate_performance_9k_ids()

        # Print summary and recommendations
        print_summary_and_recommendations()

    except KeyboardInterrupt:
        print("\n\nTests interrupted by user")
    except Exception as e:
        print(f"\n\nFatal error: {e}")
        import traceback

        traceback.print_exc()


def print_summary_and_recommendations():
    """Print final summary and recommendations."""
    print(f"\n{'='*80}")
    print(f"SUMMARY AND RECOMMENDATIONS")
    print(f"{'='*80}\n")

    # Check which approaches succeeded
    batch_success = RESULTS.get("batch_efetch", {}).get("success", False)
    gene2refseq_success = (
        RESULTS.get("bulk_files", {}).get("gene2refseq", {}).get("success", False)
    )
    gene2accession_success = (
        RESULTS.get("bulk_files", {}).get("gene2accession", {}).get("success", False)
    )
    epost_success = RESULTS.get("epost_efetch", {}).get("success", False)

    print("✓ = Successful, ✗ = Failed\n")
    print(f"  {'✓' if batch_success else '✗'} Batch efetch (10 IDs/request)")
    print(f"  {'✓' if gene2refseq_success else '✗'} gene2refseq.gz bulk catalog")
    print(f"  {'✓' if gene2accession_success else '✗'} gene2accession.gz bulk catalog")
    print(f"  {'✓' if epost_success else '✗'} epost + efetch strategy")

    print("\n" + "=" * 80)
    print("RECOMMENDED APPROACH:")
    print("=" * 80 + "\n")

    # Determine best approach based on results
    if gene2refseq_success or gene2accession_success:
        best_bulk = "gene2refseq" if gene2refseq_success else "gene2accession"
        coverage = RESULTS["bulk_files"][best_bulk].get("coverage", 0)

        print(f"PRIMARY: Use {best_bulk}.gz bulk catalog")
        print(f"  - Coverage: {coverage:.1f}% of test IDs")
        print(f"  - Speed: ~instant lookups after one-time download")
        print(f"  - Implementation: Download once, build dict/SQLite index")
        print(f"  ")
        print(f"FALLBACK: Batch efetch for remaining ~{100-coverage:.1f}% of IDs")

        if batch_success:
            time_per = RESULTS["batch_efetch"].get("time_per_record", 1.0)
            print(f"  - Speed: {time_per:.3f}s per ID (vs 1.0s for individual)")
            print(f"  - Batch size: 10-50 IDs per request")

        print("\nIMPLEMENTATION PLAN:")
        print("  1. Download + index gene2refseq.gz at pipeline setup")
        print("  2. First pass: bulk lookup all IDs (instant)")
        print("  3. Second pass: batch efetch for unresolved IDs")
        print("  4. Expected total time: < 5 minutes for 9,000 IDs")

    elif batch_success:
        time_per = RESULTS["batch_efetch"].get("time_per_record", 1.0)
        speedup = 1.0 / time_per

        print(f"Recommended: Batch efetch")
        print(f"  - Speed: {time_per:.3f}s per ID ({speedup:.1f}x speedup)")
        print(f"  - Batch size: 10-50 IDs per request")
        print(f"  - Expected time for 9,000 IDs: ~{9000*time_per/60:.1f} minutes")

    elif epost_success:
        print(f"Recommended: epost + efetch")
        print(f"  - Can handle larger batches (100+ IDs)")
        print(f"  - More complex implementation")

    else:
        print("WARNING: All alternative approaches failed!")
        print("  - May need to stick with individual efetch")
        print("  - Consider running tests with NCBI API key for better results")


if __name__ == "__main__":
    main()
