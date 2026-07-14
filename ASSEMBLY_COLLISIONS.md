# Assembly / provenance collisions — working doc

_Date: 2026-07-14 · Source: `results/resolved_ids.tsv` (16,024 resolved rows, run 2026-07-13)_

For species whose transcripts come from several provenances (resolver streams),
the pipeline can assign rows a shared or un-fetchable `assembly_accession`,
causing wrong-genome extraction or guaranteed failures. This doc inventories
those collisions so we can work through them.

## How extraction keys the genome

`extract_sequences` maps each row to disk by
`resources/cache/<assembly_accession>/genome.fasta`
([extract_sequences.py:129](workflow/scripts/extract_sequences.py#L129)). So the
unit of collision is the **`assembly_accession` value**, not the transcript.

**Merge is clean per-transcript:** every `transcript_id` is unique in the merged
output; 0 transcripts carry conflicting assemblies. The problems below are all at
the accession→cache level.

---

## 🔴 Class 1 — one `assembly_accession` shared across species (silent wrong-FASTA)

One cache dir feeds multiple, genuinely different genomes.

| `assembly_accession` | species | rows | verdict |
|---|---|---|---|
| **`Phytozome`** | **9** (amborella, chlamydomonas, citrus, oryza, physcomitrella, ricinus, sorghum, vitis, zea) | **503** | **CRITICAL** |
| `GCF_000001765.3` | 2 name-variants of *D. pseudoobscura* | 67 | benign (same GCF) |
| `GCF_000146045.2` | 2 name-variants of *S. cerevisiae* | 65 | benign (same GCF) |
| `GCF_000219495.3` | 2 name-variants of *M. truncatula* | 431 | benign (same GCF) |
| `GCF_000754195.1` | 2 name-variants of *D. simulans* | 61 | benign (same GCF) |

**`Phytozome` is a latent corruption bomb.** All 9 phytozome species collapse to
`resources/cache/Phytozome/genome.fasta`. Today it fails *safe* (no phytozome
FASTA is ever downloaded → `assembly_not_cached`), but the instant any one genome
lands there, all 9 species get sliced against it → wrong sequences, silently.

- **Status: FIXED** by commit `d546a0e` — the phytozome resolver now emits
  `assembly_accession = phytozome_<species>`, and `download_phytozome_fasta`
  caches per-species genomes. Re-run required for it to take effect in
  `resolved_ids.tsv`.
- **Deploy caveat:** delete any stale `resources/cache/Phytozome/` before/after
  the switch so extraction can't read an orphaned wrong genome.

The 4 benign rows are the *same* GCF under two organism spellings
(`medicago truncatula` vs `…(barrel medic)`). Harmless for extraction; a
naming-hygiene smell (organism strings not normalized before grouping).

---

## 🔴 Class 2 — same build, un-fetchable accession from one provenance (guaranteed fail)

A species is resolved by several streams; one emits the real downloadable GCF,
another emits a **bare build name or sequence-level accession** no downloader can
fetch — even though a sibling stream already knows the GCF for that exact build.

| provenance | bad `assembly_accession` | species | rows | real GCF (already found by NCBI stream) |
|---|---|---|---|---|
| ensembl | `ARS-UCD2.0` | *B. taurus* | 53 | `GCF_002263795.3` |
| flybase | `BDGP6` | *D. melanogaster* | 33 | `GCF_000001215.4` |
| wormbase | `WBcel235` | *C. elegans* | 17 | `GCF_000002985.6` |
| ncbi | `NT_033777.3`, `NT_033778.4`, `NT_033779.5`, `NT_037436.4`, `NT_479536.1` | *D. melanogaster* / *D. simulans* | 227 | `GCF_000001215.4` (dmel) |

These are exactly the report's ensembl-53 / flybase-33 `assembly_not_cached`
failures (~330 rows total, guaranteed recoverable).

**Root causes / fix locations:**
- **ensembl** — `resolve_ensembl_assembly_accessions.py` `ASSEMBLY_NAME_MAPPING`
  has no `ARS-UCD2.0 → GCF` entry, so bos taurus (and any build not in the static
  map) leaks the bare name. Add the missing build→GCF entries.
- **flybase / wormbase** — map *most* rows to the GCF but leak a subset as the
  bare build name; the name→GCF conversion is incomplete/inconsistent within the
  resolver. Make it total.
- **ncbi `NT_` leakage** — sequence-level RefSeq accessions (chromosome arms)
  reach extraction as `assembly_accession` instead of being mapped to the parent
  assembly GCF. Gap in the NCBI chromosome/assembly-accession resolution.

---

## 🟡 Class 3 — same species split across multiple genome builds (consistency, not a bug)

~30 species mix 2–7 builds. Each row extracts from its *own* build correctly
(per-row `acc↔coords` are consistent), so there is **no wrong-FASTA** — but the
output for a given species blends builds. This is inherent to multi-source input;
the decision is whether the downstream classifier needs one canonical build per
species.

Examples: *D. rerio* GRCz10/GRCz11/GRCz12tu · *G. max* v2.0/v2.1/v4.0 ·
*P. abelii* 4 builds · *G. gorilla* 5 builds · *S. lycopersicum* SL3.0/SLM_r2.1/
GCA_000188115.2 · *Z. mays* B73_RefGen_v4/NAM-5.0/Phytozome.

---

## zea_mays — provenance mismatch (from the phytozome-FASTA work)

`zea_mays`'s phytozome config `gtf:` is an **ensembl_plants AGPv4 GFF3**
(`resources/ensembl_plants/Zea_maysb73v4.AGPv4.gff3.gz`), but
`download_phytozome_fasta` would pull JGI **genome_id 833**'s assembly. Coords
(AGPv4, ensembl chrom names) won't match a JGI FASTA → `chrom_not_found`, not
corruption. Only 2 rows. Fix: point zea at a real JGI phytozome annotation for
833, or route it through the ensembl_plants stream instead of phytozome.

---

## Suggested order of work

1. **Class 1 / Phytozome** — already fixed in code (`d546a0e`); just guard the
   stale `cache/Phytozome/` dir on the next run. ✅ mostly done
2. **Class 2** — ~330 recoverable rows, pure mapping-completeness fixes in the
   ensembl / flybase / wormbase / ncbi assembly resolvers. Highest ROI next step.
3. **Class 3 + zea_mays** — judgment calls (canonical build per species; zea
   config routing).

## Data provenance

Reproduce the tables from `results/resolved_ids.tsv`:
- Class 1: group by `assembly_accession`, count distinct normalized `organism`.
- Class 2: rows whose `assembly_accession` is not `GC[FA]_…` / `url_…` / `phytozome_…`.
- Class 3: group by normalized `organism`, count distinct `assembly_accession`.
