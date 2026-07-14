# Phytozome genome-FASTA extraction — design

_Date: 2026-07-14 · Branch target: `main` · Related audit: [PIPELINE_AUDIT.md](../../../PIPELINE_AUDIT.md) Priority 0/1_

## Problem

503 phytozome transcripts resolve to coordinates but fail extraction with
`assembly_not_cached` — the single largest extraction-failure bucket. The
phytozome stream downloads GFF3 **annotations** (for coordinate resolution) but
never a genome **FASTA**, so `extract_sequences` has nothing to slice.

Two concrete blockers today:

1. **No genome FASTA is fetched for phytozome.** JGI concerns (bearer auth,
   PURGED tape restore) live only in the GFF3 download path.
2. **Cache-key collision.** [resolve_phytozome_gtf.py](../../../workflow/scripts/resolve_phytozome_gtf.py)
   sets `assembly_accession = "Phytozome"` for every row. Extraction maps a row
   to disk via `CACHE_DIR / assembly_accession / genome.fasta`
   ([extract_sequences.py:129](../../../workflow/scripts/extract_sequences.py#L129)),
   so all species would collide on one cache dir.

## Goal

Extract flanking sequences for phytozome-resolved transcripts by fetching the
genome FASTA of the **same Phytozome version** the coordinates were resolved
against, and caching it where extraction already looks.

## Version-matching authority (decided)

The genome FASTA MUST be the mate of the exact GFF3 used for resolution.
Authority = **the resolved annotation's JGI filename**, not config
`phytozome_version` (which stays a human note) and not input-ID parsing.

Match by filename stem: strip `.gene_exons.gff3.gz` / `.gene.gff3.gz` /
`.gff3.gz` from the resolved annotation `file_name`
(e.g. `Csinensis_154_v1.1.gene.gff3.gz` → `Csinensis_154_v1.1`), then pick the
assembly `.fa.gz` whose name shares that stem, deprioritizing masked assemblies
(`masked`/`softmasked`/`hardmasked`). If no stem match, fall back to the current
best-assembly heuristic (versioned, unmasked) with a logged warning.

## Approach (chosen: dedicated FASTA rule, JGI concerns contained)

Rejected alternatives:
- Route the FASTA through the generic `download_assembly` manifest — forces JGI
  bearer auth + restore into the unauthenticated generic downloader
  ([download_manifest_utils.py:29-31](../../../workflow/scripts/download_manifest_utils.py#L29-L31)
  also gates the manifest to NCBI-accession-or-`fasta_url`, excluding phytozome).
- One combined rule producing GFF3 + FASTA — downloads multi-GB assemblies (and
  triggers tape restores) for every configured species, even those with zero
  resolved rows.

### 1. Reusable JGI lookup — [jgi_phytozome_lookup.py](../../../workflow/scripts/jgi_phytozome_lookup.py)

- Add version-matched assembly selection keyed off a resolved annotation's
  `file_name` (stem match above; masked deprioritized; graceful fallback).
- New `resolve_pair(genome_id, token, prefer_name=None, per_page=...)` returning
  `{"annotation": describe(...), "sequence": describe(...)}` where `sequence` is
  the version-matched assembly. Thin `resolve_sequence(...)` wrapper.
- `resolve_annotation()` unchanged — the GFF3 rule is untouched.

### 2. Per-species cache key — [resolve_phytozome_gtf.py](../../../workflow/scripts/resolve_phytozome_gtf.py)

- One-line change: `assembly_accession = f"phytozome_{species}"` (species = config
  key, already the resolver's per-row species) instead of the shared
  `"Phytozome"`. Extraction then finds `resources/cache/phytozome_<species>/genome.fasta`
  with **no change to `extract_sequences.py`**.

### 3. New download rule + script — `download_phytozome_fasta`

Mirrors `download_phytozome_gtf`:
- Load JGI token; look up `genome_id` + `portal_file_name` from
  `phytozome_gtf_sources` by species (the `{species}` wildcard).
- `resolve_pair(genome_id, token, prefer_name=portal_file_name)` → assembly FASTA.
- Download `.fa.gz` with the bearer token, gunzip → `genome.fasta`,
  `samtools faidx` → `genome.fasta.fai`.
- Output `resources/cache/phytozome_{species}/{genome.fasta, genome.fasta.fai, .download_done}`.

**Fault tolerance — follow `download_assembly`, not the GFF3 rule.** Always exit
0 and write a `.download_done` sentinel (`ok` / `failed: <reason>`); on PURGED,
fire `request_restore` and record `failed: restore requested`. A cold-storage or
missing FASTA then surfaces as extraction `assembly_not_cached` (unchanged
behavior) instead of blocking the whole pipeline (the GFF3 rule blocks because
it is upstream of resolution; the FASTA is only needed for extraction).

gunzip + `samtools faidx` is a ~6-line duplication of logic in
[download_assembly.py](../../../workflow/scripts/download_assembly.py); inline it
in the new script with a `ponytail:` comment rather than extract a shared module
for two call sites.

### 4. Fan-out wiring — new rule file, mirrors [download_assemblies.smk](../../../workflow/rules/download_assemblies.smk)

- `get_phytozome_species(wc)`: gate on the existing `prepare_accession_list`
  checkpoint (guarantees the resolved table is final), then read distinct
  `db_source == phytozome` species from **`results/resolved_ids.tsv`** (merged;
  post-merge so coordinate-less rows dropped by merge are excluded; provably the
  same phytozome species set as the file extraction consumes).
- Add `expand("resources/cache/phytozome_{species}/.download_done", species=get_phytozome_species(wc))`
  as an input to `download_assemblies_done`, so the existing `.assemblies_ready`
  sentinel waits for phytozome FASTAs. Extraction is untouched.
- Only species with resolved rows get a FASTA downloaded (not all ~11 configured).

### 5. Testing — [tests/test_jgi_phytozome_lookup.py](../../../tests/test_jgi_phytozome_lookup.py)

Extend with the version-matched assembly selection (pure, offline):
- exact stem match picks the correct-version `.fa.gz`;
- masked assemblies deprioritized when an unmasked mate of the same stem exists;
- no stem match → documented fallback to best versioned/unmasked assembly.

## Data flow

```
resolve_phytozome_gtf  ──sets assembly_accession = phytozome_<species>──┐
                                                                        │
merge_resolved → resolved_ids.tsv → resolve_ncbi_chromosome_accessions  │
        │                                                               │
        └─ get_phytozome_species (gated on prepare_accession_list) ─────┤
                                                                        ▼
download_phytozome_fasta (×species) ─ resources/cache/phytozome_<species>/genome.fasta{,.fai}
                                                                        │
                            (input to) download_assemblies_done → .assemblies_ready
                                                                        │
                                                          extract_sequences (unchanged)
```

## Out of scope / deferred follow-ups

- **Rename `ncbi_chromosome_resolved.tsv`** to something explanatory (e.g.
  `resolved_chrom_patched.tsv`): it is the full merged resolved table with NCBI
  chromosome names patched, not an NCBI-only file. Touches
  `resolve_ncbi_chromosome_accessions.smk`, `download_assemblies.smk`,
  `extract_sequences.smk`. Tracked as a later refactor (tag with a `ponytail:`
  note at a reference site during implementation).
- ensembl (53) / flybase (33) / noncode (34) `assembly_not_cached` residuals —
  separate URL-gap work, not this change.
- The 71 `plant_gtf` `chrom_not_found` — separate chromosome-naming stream.

## Verification

- Unit tests above pass offline.
- `python workflow/scripts/jgi_phytozome_lookup.py --genome-id 154 --json` shows a
  version-matched `sequence` alongside the `annotation`.
- After a run: phytozome `assembly_not_cached` count drops (RESTORED species
  extract; PURGED species report `failed: restore requested` until staged);
  `grep -c '^>' results/output.fasta` rises by the recovered phytozome rows.
