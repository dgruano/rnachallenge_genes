# Design — Tool-aware plant transcript resolution (RNAPlonc / PreLnc)

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

Out of scope (deferred, need source articles): PLncPRO / CNIT maize & vitis; the exact
LGC-rice vs PreLnc-rice tie-break beyond what namespace+tool already give.

## 4. Architecture — add the tool as a resolution key

**Principle:** namespace already disambiguates *most* cases (GRMZM vs Zm00001d split maize
without the tool). The tool is required only where one namespace is shared across tools with
different sources — today that is **rice `Os…`** (PreLnc v44 vs LGC r30). So we add a
*minimal* tool-aware layer, not a rewrite.

**4a. Join the tool into the plant stream.** In the plant/phytozome resolvers (or a small
shared helper), left-join `tool_source_map.tsv` on `transcript_id` to attach `primary_tool`.
This is the only new data dependency.

**4b. Config gains an optional per-tool override.** Extend the per-species config entry with
an optional `by_tool:` map. Species with no `by_tool` behave exactly as today (backward
compatible). Example:

```yaml
oryza_sativa:
  # default (LGC etc.) stays on the existing release
  url: ".../release-30/.../Oryza_sativa.IRGSP-1.0.30.gtf.gz"
  release: "30"
  by_tool:
    PreLnc: { url: "FILL_ME__ensembl_plants_v44_oryza_gtf", release: "44" }
```

Resolver rule: `source = entry.by_tool.get(primary_tool, entry_default)`.

**4c. Multi-tool transcripts (OPEN DECISION — needs your call).** A transcript can list
several tools (e.g. `GRMZM2G703059_T01` is in both RNAPlonc and PreLnc). `dani_notes.md`
muses "maybe assume they all come from PreLnc". Proposed default precedence when a namespace
is shared: **namespace wins first** (GRMZM is AGPv3 regardless of tool), and only for a
genuinely tool-ambiguous namespace (rice `Os`) do we apply a tool precedence order
(proposed: `PreLnc > LGC`). Please confirm or override.

## 5. Config changes

- **`phytozome_gtf_sources.yaml`** — repoint citrus/potato/oryza(maize?) to v11; **add
  `brachypodium_distachyon` + `manihot_esculenta`**. `genome_id` / `portal_file_name` for
  the added/repointed species resolved via `jgi_phytozome_lookup.py` (JGI token present in
  `.env`). PURGED files fire a restore per existing CLAUDE.md workflow.
- **`plant_gtf_sources.yaml`** (or a new `ensembl_plants_v44` block) — **placeholder** v44
  entries for the PreLnc-coding namespaces: maize `Zm00001d` (AGPv4), rice `Os` (IRGSP-1.0),
  arabidopsis (TAIR10). URLs marked `FILL_ME__…` for you to drop in.

## 6. Resolver changes

- Attach `primary_tool` (§4a) in `resolve_plant_gtf.py` and/or `resolve_phytozome_gtf.py`.
- Honor `by_tool` override (§4b).
- Add missing prefixes only if a namespace routes to a resolver that lacks it (verify with a
  dry-run; `Bradi/Manes/Sobic` already classify as `plant` + phytozome `species_hint`).

## 7. Audit + rerun

- Update `PIPELINE_AUDIT.md`: supersede L6's "AGPv3/v44" guesswork with this tool-aware
  model; record the corrected RNAPlonc=Phytozome-v11 / PreLnc=v44+GreeNC mapping.
- Rerun: regenerate `classified_ids.tsv` → plant/phytozome resolvers → merge → download →
  extract → report. Because `parse_ids.py` is Stage 1, a namespace/routing change forces a
  broad rerun (noted in audit TB-3). Exact command supplied after implementation.

## 8. Explicitly NOT doing

- No GreeNC downloader (substituted by Phytozome v11).
- No PLncPRO/CNIT/vitis source wiring (needs articles; deferred).
- No change to non-plant PreLnc (human/mouse/cattle/RefSeq resolve via other paths).
