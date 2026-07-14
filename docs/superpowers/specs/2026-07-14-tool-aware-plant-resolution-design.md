# Design — Namespace-keyed plant transcript resolution (RNAPlonc / PreLnc)

_(Earlier drafts explored a tool-aware layer; §4 records why it was dropped for
namespace-keyed resolution + ID-presence verification.)_

_Date: 2026-07-14 · Branch: `fix/rnaplonc-ensembl-plants-v44` · Status: DRAFT for review_

## 1. Problem & corrected source model

The initial framing ("Ensembl Plants v44 for all RNAPlonc transcripts") was wrong on
the database. Ground truth from `dani_notes.md` + user clarification:

| Tool | Coding source | Non-coding source |
|---|---|---|
| **RNAPlonc** | Phytozome v11 (GreeNC uses Phytozome V10; V10 API had no annotation files, so **v11** is the substitute) | Phytozome v11 |
| **PreLnc** | **Ensembl Plants v44** | GreeNC (≈ Phytozome) |

Two consequences:
1. **GreeNC is not a new resolver.** GreeNC is built on Phytozome, and `dani_notes.md`
   already substitutes it with the existing **Phytozome v11** path. No GreeNC download code.
2. **Resolution must key on `(species, tool)`, not species alone.** The same species is
   served by different assemblies depending on the tool the transcript came from. Maize is
   the clearest case: `Zm00001d…` (PreLnc coding) → Ensembl Plants v44 / AGPv4, but
   `GRMZM…` (RNAPlonc) → Phytozome v11 / AGPv3. Rice `Os…t…_0N` is shared between PreLnc
   (v44) and LGC (release-26/30) under **one namespace**, so there the tool is the *only*
   discriminator.

## 2. What the pipeline does today (confirmed in code)

- `parse_ids.py` → every plant ID gets `db_source = plant`, plus `species_hint` and a
  `source_hint` (`phytozome` / `ensembl_plants` / `tair`…). `source_hint` is descriptive
  only; nothing routes on it.
- `resolve_plant_gtf.py` and `resolve_phytozome_gtf.py` both consume `db_source == plant`
  rows and resolve **by species only** (`resolve_phytozome_gtf` claims any row whose
  `species_hint` is a key in `phytozome_gtf_sources.yaml`).
- The tool is known per transcript in `results/tool_source_map.tsv`
  (`transcript_id → tools`, `primary_tool`) but is **never consulted during resolution**.
- Config today is species-keyed and version-mixed: `phytozome_gtf_sources.yaml` has
  sorghum at v11 but citrus/potato/oryza at v10 and maize pointed at an Ensembl AGPv4 file;
  `brachypodium_distachyon` and `manihot_esculenta` are **absent** (→ RNAPlonc `Bradi`/`Manes.`
  currently unresolved, matching `article_notes.md` `missing_coordinates`).

## 3. Target source matrix (this branch)

| Species | Namespace | Tool(s) | Source | Assembly | Action |
|---|---|---|---|---|---|
| maize | `GRMZM…` | RNAPlonc, PreLnc-nc | Phytozome v11 | RefGen_v3 (AGPv3) | add/repoint Phyto v11 maize |
| maize | `Zm00001d…` | PreLnc coding, PLncPRO | Ensembl Plants v44 | AGPv4 | v44 entry |
| rice | `Os…t…_0N` | PreLnc coding | Ensembl Plants v44 | IRGSP-1.0 | v44 entry |
| rice | `Os…t…_0N` | LGC | Ensembl release-26/30 | IRGSP-1.0 | already configured (LGC path) |
| arabidopsis | `AT…`,`ATCG…` | PreLnc | Ensembl Plants v44 | TAIR10 | v44 entry |
| sorghum | `Sobic.…` | RNAPlonc | Phytozome v11 | v3.1 | ✅ already v11 |
| brachypodium | `Bradi…` | RNAPlonc | Phytozome v11 | Bd21 | **ADD (JGI lookup)** |
| manihot | `Manes.…` | RNAPlonc | Phytozome v11 | v? | **ADD (JGI lookup)** |
| citrus | `orange…` | RNAPlonc | Phytozome v11 | v1.1 | repoint v10→v11 |
| potato | `PGSC…` | RNAPlonc | Phytozome v11 | v? | repoint v10→v11 |

Out of scope (deferred, need source articles): PLncPRO / CNIT maize & vitis. Rice `Os…` needs
no tie-break — PreLnc(v44) and LGC(r30) share the IRGSP-1.0 assembly, so coordinates agree.

## 4. Architecture — namespace-keyed resolution, verified by ID presence

**Decision (2026-07-14, user-approved):** resolve by **namespace → source**; the tool is
*not* a resolution key. The `by_tool` layer from the earlier draft is **dropped**.

**Why the tool is redundant — proven empirically, not asserted.** The worry was that one
species is served by different assemblies across tools. Tested by ID presence against an
on-disk AGPv4 maize annotation:

- **0 of 100** input `GRMZM…` cores are findable in AGPv4 (`Zea_maysb73v4.AGPv4.gff3.gz`),
  even with generous substring matching. `Zm00001d…` is the native AGPv4 namespace (2.8M
  lines). → RNAPlonc's `GRMZM` and PreLnc-coding's `Zm00001d` are **disjoint namespaces on
  different assemblies**. The tool split *is* the namespace split.
- Rice `Os…t` shared by PreLnc(v44) + LGC(r30): both annotate the **same assembly**
  (IRGSP-1.0), so genomic coordinates agree; the tool only selects an annotation *release*.

No plant namespace in this input maps to two assemblies. So namespace suffices.

**The correctness guarantee is verification, not trust.** We do not rely on the
namespace→assembly table being right. Each plant resolver computes a **per-species match
rate** and **fails loudly** (non-zero exit) when a species that has input IDs matches
≈0 against its configured source — the exact signature of a mis-pointed assembly (e.g.
`GRMZM` accidentally aimed at AGPv4). A `min_match_rate` threshold (default e.g. 0.05, and
a `--no-strict` escape hatch) gates the assertion. This turns a silent zero-coordinate bug
into an immediate, named failure.

## 4bis. Differentiating assembly versions (the note)

Same species ≠ same assembly. Tell them apart, in priority order:

1. **Read the header's inline assembly tag when present.** PreLnc-coding maize headers carry
   it literally: `>Zm00001d014535_T001 cdna chromosome:B73_RefGen_v4:5:…` → **AGPv4**.
2. **Else the ID prefix/format encodes the assembly generation.** Maize is the canonical case:

   | ID shape | Assembly | Ensembl/Phytozome home |
   |---|---|---|
   | `GRMZM2G…_T0N` / `AC……_FG…` | B73 **RefGen_v3** (AGPv3) | Phytozome v11 / GreeNC |
   | `Zm00001d……_T00N` | B73 **RefGen_v4** (AGPv4) | Ensembl Plants **v44** |
   | `Zm00001e……` | B73 **NAM 5.0** | Ensembl Plants ≥ r49 |

   Rice: `Os…t…_0N` (RAP/IRGSP-1.0) vs `LOC_Os…` (MSU7) — same IRGSP-1.0 genome, different
   annotation namespace. Arabidopsis `AT…G…` → TAIR10.
3. **Always confirm by ID presence** before trusting a source:
   `zgrep -c -Ff <sample_ids> <candidate>.gff3.gz`. A near-zero count means wrong assembly.
   This is the same check the resolver automates (§4).

**Consequence for config keying:** because `GRMZM` and `Zm00001d` share `species_hint =
zea_mays` but need different assemblies, maize cannot be a single species entry. It is split
by namespace: `GRMZM` → Phytozome v11 maize (AGPv3); `Zm00001d` → Ensembl Plants v44 (AGPv4).
See §5 for how the split is keyed.

## 5. Config changes

**RNAPlonc → Phytozome v11** (`phytozome_gtf_sources.yaml`):
- **Add** `brachypodium_distachyon` (`Bradi…`) + `manihot_esculenta` (`Manes.…`) — currently
  unresolved. `genome_id` / `portal_file_name` resolved via `jgi_phytozome_lookup.py` (JGI
  token in `.env`); PURGED files fire a restore per CLAUDE.md.
- **Maize `GRMZM` split:** add a Phytozome v11 maize (AGPv3) entry so `GRMZM` stops being
  aimed at the AGPv4 file. Because the resolver keys on `species_hint`, the maize split is
  enforced in `parse_ids.py`: give `Zm00001d…` a distinct `species_hint`
  (`zea_mays` stays for `GRMZM`→Phytozome; `Zm00001d`→a v44 entry, §below). Verified by the
  presence check (§4) — GRMZM must match the Phytozome file, not AGPv4.
- Citrus/potato already resolve on their current Phytozome files; leave version as-is unless
  the presence check flags a miss (do **not** churn working entries — ponytail).

**PreLnc-coding → Ensembl Plants v44** (`plant_gtf_sources.yaml`, `release: "44"` entries with
**`FILL_ME__` placeholder URLs** for you to fill):
- rice `oryza_sativa` v44 (IRGSP-1.0), arabidopsis v44 (TAIR10), maize `Zm00001d` v44 (AGPv4).
- Note: `Zm00001d` headers already carry inline coordinates (`chromosome:B73_RefGen_v4:…`);
  the v44 GTF is the coordinate source of record and lets the presence check validate them.

## 6. Resolver changes

- **Presence-verification guard** (the correctness core): in `resolve_plant_gtf.py` and
  `resolve_phytozome_gtf.py`, after resolving, compute per-species match rate =
  matched / (matched + unmatched with a configured source). If any such species falls below
  `min_match_rate`, log an ERROR listing the species+source and exit non-zero (unless
  `--no-strict`). One small unit test on the rate/threshold logic.
- **Routing:** enforce the maize `GRMZM`/`Zm00001d` `species_hint` split in `parse_ids.py`;
  confirm `Bradi/Manes/Sobic` reach the phytozome resolver via a dry-run (they classify as
  `plant` + phytozome `species_hint` already). Add prefixes only where a dry-run shows a gap.

## 7. Audit + rerun

- Update `PIPELINE_AUDIT.md`: supersede L6's AGPv3/v44 guesswork with the verified
  namespace→assembly model; record RNAPlonc=Phytozome-v11 / PreLnc=v44+GreeNC and the
  presence-check guard.
- Rerun: regenerate `classified_ids.tsv` → plant/phytozome resolvers → merge → download →
  extract → report. Because `parse_ids.py` is Stage 1, a namespace/routing change forces a
  broad rerun (noted in audit TB-3). Exact command supplied after implementation.

## 8. Explicitly NOT doing

- No GreeNC downloader (substituted by Phytozome v11).
- No PLncPRO/CNIT/vitis source wiring (needs articles; deferred).
- No change to non-plant PreLnc (human/mouse/cattle/RefSeq resolve via other paths).
