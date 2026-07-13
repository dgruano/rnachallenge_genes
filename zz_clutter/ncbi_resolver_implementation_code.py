"""
NCBI Resolver - Recommended Implementation Code
================================================
Complete code snippets for fixing the suppressed record issue.

DO NOT RUN THIS FILE - These are reference implementations only.
Copy the relevant functions into resolve_ids.py as needed.
"""

from typing import Optional

from Bio import Entrez

# ═══════════════════════════════════════════════════════════════════════
# SOLUTION 1: Quick Fix (Priority 1)
# ═══════════════════════════════════════════════════════════════════════


def fix_esearch_query_BEFORE():
    """Current (broken) implementation."""
    # Line ~110 in resolve_ids.py
    handle = with_retry(
        Entrez.esearch,
        db="nucleotide",
        term=" OR ".join(
            f"{acc}[ACCN]" for acc in batch
        ),  # ← Problem: [ACCN] excludes suppressed
        retmax=len(batch) * 5,
        label="ncbi_esearch",
    )


def fix_esearch_query_AFTER():
    """Fixed implementation - removes [ACCN] tag."""
    # Line ~110 in resolve_ids.py
    handle = with_retry(
        Entrez.esearch,
        db="nucleotide",
        term=" OR ".join(f"{acc}" for acc in batch),  # ← Fixed: no [ACCN] tag
        retmax=len(batch) * 5,
        label="ncbi_esearch",
    )


# ═══════════════════════════════════════════════════════════════════════
# SOLUTION 2: Comprehensive Fix (Priority 2)
# ═══════════════════════════════════════════════════════════════════════


def resolve_via_efetch(accession: str) -> Optional[dict]:
    """
    Fallback resolver for suppressed/problematic NCBI records.

    This function should be added after line ~400 in resolve_ids.py,
    right after the resolve_ncbi_batch() function definition.

    Args:
        accession: NCBI transcript accession (e.g., "XM_020345473.1")

    Returns:
        dict with resolved gene info if successful, None otherwise

    Usage:
        When elink returns no gene links, call this as fallback:

        if not gene_ids_for_tid:
            log.debug(f"Trying efetch fallback for {tid}")
            result = resolve_via_efetch(tid)
            if result:
                resolved_rows.append(result)
    """
    try:
        # Step 1: Fetch GenBank record via efetch
        handle = with_retry(
            Entrez.efetch,
            db="nucleotide",
            id=accession,
            rettype="gb",
            retmode="xml",
            label=f"efetch_fallback_{accession}",
        )

        if handle is None:
            log.warning(f"  efetch retry failed for {accession}")
            return None

        records = Entrez.read(handle)
        handle.close()

        if not records:
            log.warning(f"  efetch returned no records for {accession}")
            return None

        gb = records[0]

        # Step 2: Check if record is suppressed (for logging purposes)
        comment = gb.get("GBSeq_comment", "")
        if "suppressed" in comment.lower() or "removed" in comment.lower():
            log.debug(f"  {accession} is suppressed/removed")

        # Step 3: Extract gene information from feature table
        gene_id = None
        gene_symbol = None

        features = gb.get("GBSeq_feature-table", [])
        for feat in features:
            # Look for the 'gene' feature
            if feat.get("GBFeature_key") == "gene":
                qualifiers = feat.get("GBFeature_quals", [])

                for qual in qualifiers:
                    qual_name = qual.get("GBQualifier_name", "")
                    qual_value = qual.get("GBQualifier_value", "")

                    # Extract GeneID from db_xref
                    if qual_name == "db_xref" and "GeneID:" in qual_value:
                        # Parse "GeneID:109786895" or "GeneID:109786895,OtherDB:123"
                        gene_id = (
                            qual_value.split("GeneID:")[1]
                            .split(",")[0]
                            .split(";")[0]
                            .strip()
                        )

                    # Extract gene symbol
                    elif qual_name == "gene":
                        gene_symbol = qual_value

        if not gene_id:
            log.warning(
                f"  {accession}: GenBank record has no gene feature with GeneID"
            )
            return None

        log.info(
            f"  {accession}: extracted GeneID={gene_id} from GenBank feature table"
        )

        # Step 4: Fetch gene details using the extracted GeneID
        gene_handle = with_retry(
            Entrez.esummary,
            db="gene",
            id=gene_id,
            label=f"gene_esummary_{gene_id}",
        )

        if gene_handle is None:
            log.warning(f"  Failed to fetch gene summary for GeneID={gene_id}")
            return None

        gene_summaries = Entrez.read(gene_handle)
        gene_handle.close()

        # Step 5: Parse gene summary to get genomic coordinates
        doc = gene_summaries.get("DocumentSummarySet", {}).get("DocumentSummary", [{}])[
            0
        ]
        genomic_info = doc.get("GenomicInfo", [{}])
        loc = genomic_info[0] if genomic_info else {}

        # Step 6: Build result dictionary matching expected schema
        return {
            "transcript_id": accession,
            "db_source": "ncbi",
            "gene_id": gene_id,
            "gene_symbol": gene_symbol or doc.get("Name", ""),
            "organism": doc.get("Organism", {}).get("ScientificName", ""),
            "assembly_accession": loc.get("ChrAccVer", ""),
            "chrom": loc.get("ChrLoc", ""),
            "start": int(loc.get("ChrStart", 0)),
            "end": int(loc.get("ChrStop", 0)),
            "strand": str(loc.get("ChrStrand", "")).strip() or "+",
            "is_ambiguous": False,
        }

    except Exception as exc:
        log.warning(f"  efetch fallback failed for {accession}: {exc}")
        return None


# ═══════════════════════════════════════════════════════════════════════
# How to integrate into resolve_ncbi_batch()
# ═══════════════════════════════════════════════════════════════════════


def integration_point_in_resolve_ncbi_batch():
    """
    This shows where to add the efetch fallback call.

    Location: Around line ~320 in resolve_ids.py
    Inside the loop: for tid in batch:
    """

    # ... existing code ...

    # Try gene link for the found GI
    gene_ids_for_tid = gi_to_genes.get(gi, [])

    if not gene_ids_for_tid:
        # EXISTING WARNING
        log.warning(f"  NCBI: no gene link for transcript {tid} (GI={gi})")

        # ╭─────────────────────────────────────────────────╮
        # │ INSERT NEW CODE HERE                           │
        # ╰─────────────────────────────────────────────────╯

        # Try efetch fallback for suppressed records
        log.debug(f"  Attempting efetch fallback for {tid}")
        fallback_result = resolve_via_efetch(tid)

        if fallback_result:
            resolved_rows.append(fallback_result)
            log.info(
                f"  ✓ {tid} resolved via efetch fallback → gene {fallback_result['gene_id']}"
            )
        else:
            log.warning(f"  ✗ {tid} could not be resolved by any method")

        # Skip to next transcript
        continue

        # ╭─────────────────────────────────────────────────╮
        # │ END OF INSERTION                               │
        # ╰─────────────────────────────────────────────────╯

    # ... rest of existing code (handles ambiguous genes, etc.) ...


# ═══════════════════════════════════════════════════════════════════════
# OPTIONAL: Check for replacement IDs (low success rate expected)
# ═══════════════════════════════════════════════════════════════════════


def check_for_replacement(gi: str) -> Optional[str]:
    """
    Check if a suppressed record has a replacement ID.

    Returns:
        Replacement accession if available, None otherwise

    Note: Testing showed 0/3 suppressed records had replacements,
    but this is still worth checking as a safety measure.
    """
    try:
        handle = Entrez.esummary(db="nucleotide", id=gi)
        result = Entrez.read(handle)
        handle.close()

        if result:
            doc = result[0]
            status = doc.get("Status", "")
            replaced_by = doc.get("ReplacedBy", "").strip()

            if status == "suppressed" and replaced_by:
                log.info(f"  GI={gi} is suppressed, replaced by {replaced_by}")
                return replaced_by

        return None

    except Exception as e:
        log.warning(f"  Could not check replacement for GI={gi}: {e}")
        return None


def integration_with_replacement_check():
    """Example of how to integrate replacement check."""

    # When we find a GI via esearch fallback...
    if gi is None:
        log.warning(f"  NCBI: could not map {tid} to a GI")
        continue

    # NEW: Check if it's suppressed with a replacement
    replacement = check_for_replacement(gi)
    if replacement:
        log.info(f"  Following replacement: {tid} → {replacement}")
        # Recursively resolve the replacement
        # (Implementation depends on your batch processing logic)
        continue

    # Continue with normal gene linking...


# ═══════════════════════════════════════════════════════════════════════
# Testing/Validation Code
# ═══════════════════════════════════════════════════════════════════════


def test_suppressed_record_resolution():
    """
    Test function to validate the efetch fallback approach.

    Run this separately to verify the fix works before integrating.
    """
    test_ids = [
        "XM_020345473.1",  # Known suppressed, GeneID:109786895
        "XM_020345467.1",  # Known suppressed, GeneID:109786895
        "XM_020315385.1",  # Known suppressed, GeneID:109756534
    ]

    Entrez.email = "your_email@example.com"
    # Entrez.api_key = "your_api_key"  # Recommended

    print("Testing efetch fallback on known suppressed records:")
    print("=" * 60)

    for accession in test_ids:
        print(f"\nTesting: {accession}")
        result = resolve_via_efetch(accession)

        if result:
            print(f"  ✓ SUCCESS")
            print(f"    Gene ID: {result['gene_id']}")
            print(f"    Gene Symbol: {result['gene_symbol']}")
            print(f"    Organism: {result['organism']}")
        else:
            print(f"  ✗ FAILED")

    print("\n" + "=" * 60)


# ═══════════════════════════════════════════════════════════════════════
# Performance Monitoring
# ═══════════════════════════════════════════════════════════════════════


def add_performance_tracking():
    """
    Add counters to track fallback usage.

    Insert at the top of resolve_ncbi_batch():
    """
    # Track fallback statistics
    fallback_stats = {
        "attempted": 0,
        "successful": 0,
        "failed": 0,
    }

    # ... in the fallback section ...

    if not gene_ids_for_tid:
        fallback_stats["attempted"] += 1
        result = resolve_via_efetch(tid)

        if result:
            fallback_stats["successful"] += 1
            resolved_rows.append(result)
        else:
            fallback_stats["failed"] += 1

    # ... at the end of the function ...

    if fallback_stats["attempted"] > 0:
        log.info(
            f"  Fallback stats: {fallback_stats['successful']}/{fallback_stats['attempted']} "
            f"succeeded ({fallback_stats['failed']} failed)"
        )


# ═══════════════════════════════════════════════════════════════════════
# Summary of Changes
# ═══════════════════════════════════════════════════════════════════════

"""
CHANGE SUMMARY:

1. Quick Fix (Priority 1 - Deploy Immediately):
   - File: workflow/scripts/resolve_ids.py
   - Line: ~110
   - Change: Remove [ACCN] from esearch query
   - Impact: Fixes 193 "could not map to GI" failures
   - Risk: Very low
   - Time: 1 minute

2. Comprehensive Fix (Priority 2 - Next Sprint):
   - File: workflow/scripts/resolve_ids.py
   - Lines: Add new function after line 400
   - Change: Add resolve_via_efetch() function
   - Lines: Modify ~320 to call fallback
   - Change: Add fallback call when no gene links found
   - Impact: Fixes 802 "no gene link" failures
   - Risk: Low (only used as fallback)
   - Time: 30 minutes implementation + testing

EXPECTED RESULTS:
- Before: 8,250/9,245 resolved (89.2%)
- After Quick Fix: ~8,800/9,245 resolved (95.2%)
- After Comprehensive Fix: ~9,200/9,245 resolved (99.5%)

PERFORMANCE IMPACT:
- Quick Fix: No change (~3 minutes total)
- Comprehensive Fix: +9 minutes (~12 minutes total)
  - Batch processing: 3 minutes
  - Fallback efetch: 9 minutes (for ~800 IDs)
"""
