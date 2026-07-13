"""
NCBI Entrez API Investigation Script
=====================================
Tests different approaches to resolve transcript IDs that fail via esearch
but are accessible through the NCBI web interface.

Example problematic ID: XM_020345473.1
Web URL: https://www.ncbi.nlm.nih.gov/nuccore/XM_020345473.1?report=genbank
"""

import sys
import time

from Bio import Entrez

# Configure Entrez
Entrez.email = "test@example.com"
Entrez.api_key = None  # Will use from config if available

# Test IDs that failed in resolve_ids.log
FAILED_IDS = [
    "XM_020345473.1",
    "XM_020345467.1",
    "XM_020315385.1",
    "XM_020315948.1",
    "XM_020323553.1",
    "XM_020324435.1",
    "XM_020327859.1",
]


def print_section(title):
    """Print a formatted section header."""
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80)


def test_esearch_approaches(accession):
    """Test different esearch query formats."""
    print_section(f"Testing esearch approaches for {accession}")

    base_id = accession.split(".")[0]

    approaches = [
        (f"{accession}[ACCN]", "Full ID with [ACCN] tag"),
        (f"{base_id}[ACCN]", "Base ID with [ACCN] tag"),
        (accession, "Full ID without tag"),
        (base_id, "Base ID without tag"),
        (f"{accession}[ALL]", "Full ID with [ALL] tag"),
        (f"txid{accession}", "With txid prefix"),
    ]

    databases = ["nucleotide", "nuccore"]

    for db in databases:
        print(f"\n--- Database: {db} ---")
        for term, label in approaches:
            try:
                handle = Entrez.esearch(db=db, term=term, retmax=5)
                result = Entrez.read(handle)
                handle.close()

                id_list = result.get("IdList", [])
                count = result.get("Count", 0)

                status = "✓ FOUND" if id_list else "✗ EMPTY"
                print(
                    f"  {status:12} | {label:30} | Term: {term:30} | IDs: {id_list[:3]}"
                )

                if id_list:
                    return db, id_list[0], term  # Return first successful result

            except Exception as e:
                print(f"  ERROR      | {label:30} | {str(e)[:50]}")

            time.sleep(0.4)  # Rate limiting

    return None, None, None


def test_efetch_direct(accession):
    """Try direct efetch without esearch."""
    print_section(f"Testing efetch direct access for {accession}")

    try:
        # Try fetching directly by accession
        handle = Entrez.efetch(
            db="nucleotide", id=accession, rettype="gb", retmode="xml"
        )
        record = Entrez.read(handle)
        handle.close()

        print("✓ SUCCESS: efetch returned data")
        print(f"  Records returned: {len(record)}")

        if record:
            gb = record[0]
            print(f"\n  Record details:")
            print(f"    GBSeq_locus: {gb.get('GBSeq_locus', 'N/A')}")
            print(f"    GBSeq_length: {gb.get('GBSeq_length', 'N/A')}")
            print(f"    GBSeq_moltype: {gb.get('GBSeq_moltype', 'N/A')}")
            print(f"    GBSeq_organism: {gb.get('GBSeq_organism', 'N/A')}")
            print(f"    GBSeq_create-date: {gb.get('GBSeq_create-date', 'N/A')}")
            print(f"    GBSeq_update-date: {gb.get('GBSeq_update-date', 'N/A')}")

            # Check for gene features
            features = gb.get("GBSeq_feature-table", [])
            print(f"\n  Features found: {len(features)}")

            gene_info = []
            for feat in features:
                feat_key = feat.get("GBFeature_key", "")
                if feat_key in ["gene", "CDS", "mRNA"]:
                    qualifiers = feat.get("GBFeature_quals", [])
                    quals_dict = {}
                    for qual in qualifiers:
                        name = qual.get("GBQualifier_name", "")
                        value = qual.get("GBQualifier_value", "")
                        quals_dict[name] = value

                    if feat_key == "gene":
                        print(f"\n  GENE feature found:")
                        print(f"    gene: {quals_dict.get('gene', 'N/A')}")
                        print(f"    db_xref: {quals_dict.get('db_xref', 'N/A')}")
                        gene_info.append(quals_dict)

            # Check for comment field (may contain suppression/replacement info)
            comment = gb.get("GBSeq_comment", "")
            if comment:
                print(f"\n  Comment field:")
                print(f"    {comment[:500]}")

            # Check for secondary accessions
            secondary = gb.get("GBSeq_secondary-accessions", [])
            if secondary:
                print(f"\n  Secondary accessions: {secondary}")

            return True, gb, gene_info

    except Exception as e:
        error_msg = str(e)
        print(f"✗ FAILED: {error_msg}")

        # Check if it's a "not found" vs other error
        if "400" in error_msg or "ID list is empty" in error_msg:
            print("  → Record appears to be deleted/suppressed")
        elif "Invalid uid" in error_msg:
            print("  → Invalid accession format")
        else:
            print(f"  → Unexpected error: {error_msg[:200]}")

        return False, None, None


def test_elink_approaches(accession, gi=None):
    """Test different elink strategies."""
    print_section(f"Testing elink approaches for {accession}")

    if gi is None:
        print("  No GI provided, skipping elink tests")
        return None

    print(f"  Using GI: {gi}")

    link_commands = [
        ("neighbor", "Find related records (neighbor)"),
        ("neighbor_history", "Find related records with history"),
        (None, "Default elink (no cmd specified)"),
    ]

    for cmd, label in link_commands:
        try:
            params = {
                "dbfrom": "nucleotide",
                "db": "gene",
                "id": gi,
            }
            if cmd:
                params["cmd"] = cmd

            handle = Entrez.elink(**params)
            result = Entrez.read(handle)
            handle.close()

            print(f"\n  {label}:")
            if result:
                for link_set in result:
                    print(f"    DbFrom: {link_set.get('DbFrom', 'N/A')}")
                    print(f"    IdList: {link_set.get('IdList', [])}")

                    link_set_dbs = link_set.get("LinkSetDb", [])
                    if link_set_dbs:
                        for lsdb in link_set_dbs:
                            print(f"    LinkName: {lsdb.get('LinkName', 'N/A')}")
                            links = lsdb.get("Link", [])
                            print(f"    Gene IDs: {[l['Id'] for l in links[:5]]}")
                            if links:
                                return lsdb.get("LinkName"), [l["Id"] for l in links]
                    else:
                        print(f"    ✗ No LinkSetDb found")

            time.sleep(0.4)

        except Exception as e:
            print(f"  ✗ Error: {str(e)[:100]}")

    return None


def test_esummary_direct(accession):
    """Try esummary with accession directly."""
    print_section(f"Testing esummary with accession {accession}")

    try:
        handle = Entrez.esummary(db="nucleotide", id=accession)
        result = Entrez.read(handle)
        handle.close()

        print("✓ SUCCESS: esummary returned data")
        print(f"  Result: {result}")
        return True, result

    except Exception as e:
        print(f"✗ FAILED: {str(e)}")
        return False, None


def main():
    """Run all investigation tests."""
    print(
        """
╔══════════════════════════════════════════════════════════════════════════════╗
║  NCBI Entrez API Investigation - Transcript ID Resolution Failures         ║
╚══════════════════════════════════════════════════════════════════════════════╝

This script investigates why certain NCBI transcript IDs fail via Entrez.esearch
but are still accessible through the NCBI web interface.

Testing {count} failed IDs from resolve_ids.log
    """.format(
            count=len(FAILED_IDS)
        )
    )

    # Summary results
    results = {
        "esearch_success": [],
        "efetch_success": [],
        "esearch_failed": [],
        "efetch_failed": [],
    }

    for i, accession in enumerate(FAILED_IDS[:3], 1):  # Test first 3 to save time
        print(f"\n{'#' * 80}")
        print(f"#  TEST {i}/{min(3, len(FAILED_IDS))}: {accession}")
        print(f"{'#' * 80}")

        # Test 1: esearch approaches
        db, gi, term = test_esearch_approaches(accession)
        if gi:
            results["esearch_success"].append(accession)
            print(f"\n✓ esearch found GI: {gi} using term: {term}")
        else:
            results["esearch_failed"].append(accession)
            print(f"\n✗ esearch failed for all approaches")

        time.sleep(1)

        # Test 2: efetch direct
        efetch_ok, record, gene_info = test_efetch_direct(accession)
        if efetch_ok:
            results["efetch_success"].append(accession)
        else:
            results["efetch_failed"].append(accession)

        time.sleep(1)

        # Test 3: elink if we have a GI
        if gi:
            link_results = test_elink_approaches(accession, gi)

        time.sleep(1)

        # Test 4: esummary direct
        test_esummary_direct(accession)

        time.sleep(2)  # Be nice to NCBI servers

    # Print summary
    print_section("INVESTIGATION SUMMARY")
    print(f"\nTotal IDs tested: {min(3, len(FAILED_IDS))}")
    print(f"\nesearch successes: {len(results['esearch_success'])}")
    for acc in results["esearch_success"]:
        print(f"  ✓ {acc}")

    print(f"\nesearch failures: {len(results['esearch_failed'])}")
    for acc in results["esearch_failed"]:
        print(f"  ✗ {acc}")

    print(f"\nefetch successes: {len(results['efetch_success'])}")
    for acc in results["efetch_success"]:
        print(f"  ✓ {acc}")

    print(f"\nefetch failures: {len(results['efetch_failed'])}")
    for acc in results["efetch_failed"]:
        print(f"  ✗ {acc}")

    print_section("RECOMMENDATIONS")

    if results["efetch_success"]:
        print(
            """
✓ FINDING: efetch can retrieve records that esearch cannot find

RECOMMENDATION: Implement a fallback strategy in resolve_ids.py:
  1. Try current esearch approach first (efficient for batches)
  2. For failed IDs, use efetch to directly retrieve the GenBank record
  3. Parse gene information from the GenBank feature table
  4. Extract /gene qualifier db_xref for GeneID
  5. Use that GeneID to get gene details via esummary
        """
        )

    if not results["efetch_success"] and not results["esearch_success"]:
        print(
            """
✗ FINDING: Neither esearch nor efetch can retrieve these records

This suggests:
  - Records may be truly deleted/suppressed from Entrez
  - They may be available on web interface via cached/archived data
  - May need alternative data sources (e.g., assembly reports, RefSeq)
        """
        )

    print("\n" + "=" * 80)


if __name__ == "__main__":
    main()
