"""
scripts/gramene_resolver.py
Resolve legacy plant IDs using Gramene API
==========================================
Uses Gramene search endpoint to map legacy IDs to modern Ensembl IDs
"""

import sys
from pathlib import Path
import requests
import pandas as pd
import time
from typing import Dict, Optional, Tuple

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("gramene_resolver", snakemake.log[0])
unresolved_tsv = Path(snakemake.input[0])
resolved_output = Path(snakemake.output.resolved)
unresolved_output = Path(snakemake.output.unresolved)
cfg = snakemake.config

GRAMENE_API = "https://data.gramene.org/v69/search"
REQUEST_DELAY = 0.1  # seconds between requests
MAX_RETRIES = 3

def build_query_candidates(transcript_id: str) -> list[str]:
    """
    Build ordered Gramene query candidates from transcript ID.
    
    Examples:
        GRMZM5G842623_T01 → GRMZM5G842623
        Os03t0779300_01 → Os03t0779300
        Solyc01g005000.2.1 → Solyc01g005000.2, Solyc01g005000
    """
    candidates: list[str] = [transcript_id]

    gene_id = transcript_id
    if '_T' in gene_id:
        gene_id = gene_id.split('_T')[0]
        candidates.append(gene_id)
    elif '_' in gene_id:
        head, tail = gene_id.rsplit('_', 1)
        if tail.isdigit():
            gene_id = head
            candidates.append(gene_id)

    if '.' in gene_id:
        parts = gene_id.split('.')
        if len(parts) > 2:
            candidates.append('.'.join(parts[:-1]))
        candidates.append(parts[0])

    # de-duplicate while preserving order
    seen = set()
    uniq = []
    for candidate in candidates:
        if candidate not in seen:
            seen.add(candidate)
            uniq.append(candidate)
    return uniq


def infer_species_hint(row: pd.Series) -> str:
    if 'inferred_species' in row and pd.notna(row['inferred_species']):
        return str(row['inferred_species'])
    reason = str(row.get('reason', ''))
    if reason.startswith('biomart_no_match_'):
        return reason.replace('biomart_no_match_', '')
    return ''

def query_gramene(gene_id: str) -> Optional[Dict]:
    """
    Query Gramene search API for a gene ID.
    
    Returns:
        Dict with keys: modern_id, synonyms, taxon_id, system_name
        None if not found or error
    """
    for attempt in range(MAX_RETRIES):
        try:
            params = {"q": gene_id}
            response = requests.get(GRAMENE_API, params=params, timeout=30)
            
            if response.status_code == 200:
                data = response.json()
                num_found = data['response']['numFound']
                
                if num_found == 0:
                    return None
                
                # Take first result (highest score)
                doc = data['response']['docs'][0]
                
                return {
                    'modern_id': doc['id'],
                    'synonyms': doc.get('synonyms', []),
                    'taxon_id': doc.get('taxon_id'),
                    'system_name': doc.get('system_name'),
                    'biotype': doc.get('biotype'),
                    'num_found': num_found,
                }
            else:
                log.warning(f"Gramene API returned HTTP {response.status_code} for {gene_id}")
                
        except requests.RequestException as exc:
            log.warning(f"Gramene request failed (attempt {attempt+1}/{MAX_RETRIES}): {exc}")
            if attempt < MAX_RETRIES - 1:
                time.sleep(1 * (attempt + 1))  # Exponential backoff
        except Exception as exc:
            log.error(f"Unexpected error querying {gene_id}: {exc}")
            return None
    
    return None

def resolve_via_gramene(unresolved_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Resolve transcript IDs via Gramene API.
    
    Returns:
        (resolved_df, still_unresolved_df)
    """
    resolved_records = []
    unresolved_records = []
    
    total = len(unresolved_df)
    log.info(f"Querying Gramene for {total} IDs...")
    
    for idx, row in unresolved_df.iterrows():
        transcript_id = row['transcript_id']
        species_hint = infer_species_hint(row)

        candidates = build_query_candidates(transcript_id)
        result = None
        matched_candidate = ""
        for candidate in candidates:
            result = query_gramene(candidate)
            if result:
                matched_candidate = candidate
                break
        
        if result:
            # Map modern gene ID back to transcript coordinate
            modern_gene_id = result['modern_id']
            
            resolved_records.append({
                'transcript_id': transcript_id,
                'gene_id': modern_gene_id,
                'db_source': 'gramene',
                'species': result['system_name'],
                'original_species_hint': species_hint,
                'query_id': matched_candidate,
                'gramene_synonyms': ';'.join(result['synonyms']),
                'biotype': result['biotype'],
                'num_matches': result['num_found'],
            })
            
            if (idx + 1) % 50 == 0:
                log.info(f"Progress: {idx+1}/{total} ({len(resolved_records)} resolved)")
        else:
            unresolved_records.append({
                'transcript_id': transcript_id,
                'normalized_gene_id': candidates[-1],
                'inferred_species': species_hint,
                'reason': 'not_found_in_gramene',
            })
        
        # Rate limiting
        time.sleep(REQUEST_DELAY)
    
    log.info(f"Gramene resolution complete: {len(resolved_records)}/{total} resolved ({100*len(resolved_records)/total:.1f}%)")
    
    resolved_df = pd.DataFrame(resolved_records)
    still_unresolved_df = pd.DataFrame(unresolved_records)
    
    return resolved_df, still_unresolved_df

# ── Main ───────────────────────────────────────────────────────

log.info(f"Reading unresolved IDs from {unresolved_tsv}")
unresolved_df = pd.read_csv(unresolved_tsv, sep='\t')
log.info(f"Loaded {len(unresolved_df)} unresolved IDs")

# Show breakdown by species
if 'inferred_species' in unresolved_df.columns:
    species_counts = unresolved_df['inferred_species'].value_counts()
    log.info(f"Species breakdown:\n{species_counts}")

# Resolve via Gramene
resolved_df, still_unresolved_df = resolve_via_gramene(unresolved_df)

# Save results
resolved_output.parent.mkdir(parents=True, exist_ok=True)
unresolved_output.parent.mkdir(parents=True, exist_ok=True)

if len(resolved_df) > 0:
    resolved_df.to_csv(resolved_output, sep='\t', index=False)
    log.info(f"Saved {len(resolved_df)} resolved IDs to {resolved_output}")
else:
    # Create empty file with header
    pd.DataFrame(columns=['transcript_id', 'gene_id', 'db_source']).to_csv(
        resolved_output, sep='\t', index=False
    )
    log.warning("No IDs resolved via Gramene")

if len(still_unresolved_df) > 0:
    still_unresolved_df.to_csv(unresolved_output, sep='\t', index=False)
    log.info(f"Saved {len(still_unresolved_df)} still-unresolved IDs to {unresolved_output}")
else:
    # All resolved!
    pd.DataFrame(columns=['transcript_id', 'reason']).to_csv(
        unresolved_output, sep='\t', index=False
    )
    log.info("All IDs resolved!")

# Summary statistics
if len(resolved_df) > 0:
    log.info("\n=== Gramene Resolution Summary ===")
    log.info(f"Total queried: {len(unresolved_df)}")
    log.info(f"Resolved: {len(resolved_df)} ({100*len(resolved_df)/len(unresolved_df):.1f}%)")
    log.info(f"Still unresolved: {len(still_unresolved_df)} ({100*len(still_unresolved_df)/len(unresolved_df):.1f}%)")
    
    if 'species' in resolved_df.columns:
        species_resolved = resolved_df['species'].value_counts()
        log.info(f"\nResolved by species:\n{species_resolved}")
