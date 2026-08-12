# PROJECT STATUS LEDGER — CFDE gene-set deliverables
_Last updated: 2026-07-01. Single source of truth for what's uploaded, held, running, reserved, dropped._

## ✅ UPLOADED & VERIFIED in SDG (`/humgen/.../submissions/SDG`)
| Item | Notes | Keep? |
|---|---|---|
| scRNA cNMF — SMART-seq (human) | Allen/BICCN, clean | keep |
| scRNA cNMF — M1 (human) | Allen/BICCN, clean | keep |
| **Batch 4** (1,164 sets) | GTEx eQTL/sQTL/depleted, ClinVar, ClinGen, MitoCarta×GTEx, interop | keep |
| **Batch 5** (765 sets) | shRNA-KD regulons, labeled experimental/K562-HepG2 | keep |
| GTEx_tissue_enriched | clean | keep |
| GlyGen_glyco_genesets, glyco_genesets_GTEx_ABCD | clean (UniProt crosswalk ok) | keep |
| ENCODE_ATAC / DNase / histone_marks / rE2G | ⚠️ **absolute peak→gene, no control** | **TAKE DOWN → replaced by Batch 6 / region-level** |
| cluster_accessibility_scripts, cluster_scrna_and_4dn_pipeline | code bundles | keep (code) |

## ⏳ BUILT, HELD (not uploaded — waiting on accessibility to finish)
| Item | Sets | Why held |
|---|---|---|
| **Batch 6** — bg-corrected accessibility (gene-promoter level) | 3,567 | Superseded by region-level; **may drop** in favor of region-level |
| Batch 3 (local) — GlyGen-observed + RNA-seq (clean) + H3K9/histone×GTEx (flagged) | mixed | split: clean parts ship, histone parts → region-level |
| TF ChIP regulons (+×GTEx) | ~61.8k | held — does contrast requirement extend to TF binding? (advisor) |
| eCLIP RBP regulons (+×GTEx) | ~9.1k | held — same question for eCLIP |
| eCLIP∩shRNA-KD (not yet built) | — | held pending above |

## ✅ ACCESSIBILITY REDO — COMPLETE (advisor-approved; ready to upload, supersedes flagged items)
| Batch | Contents | Sets | MD5 |
|---|---|---|---|
| **Batch 7** | region-level accessibility (ATAC/DNase/7 histone, Up/Down) | 2,681 | 1c46aa2e9496f87d530c427cce9cf047 |
| **Batch 8** | accessibility×expression concordance (ENCODE-RNAseq+GTEx, raw+control+delta) | 2,773 | 1e22e8c8b7f753beaad959da9b950fdd |
- Kept regions (local, other projects): `accessibility_regions/regions/` = 1,732 peak BEDs + 2,746 ±region BEDs.
- Batch 6 (gene-level bgcontrast) SUPERSEDED by Batch 7 — do not upload.
- rE2G: symbol-cleaned gene-level bgcontrast only (not region-level'd; it's a link model, not raw peaks).

## scRNA (cluster) — SMART-seq ✅ingested, M1 ✅ingested, mouse-10x ✅submitted(log-confirmed); mouse-SMARTseq running (4th)

## 🔒 RESERVE (`_RESERVE_hold_for_later/`) — NOT for incubator
iPTMnet; CRISPR/CRISPRi KD (999); GENCODE cis-antisense; gnomAD constraint; HPO (5,261); gnomAD-dependent interop.

## ❌ DROPPED
Pharos/IDG — public GraphQL won't paginate (returns same 10 regardless of skip). Revisit via bulk source post-deadline.

## 📌 OPEN DECISIONS
1. Accessibility design: advisor **signed off** (background-prevalence contrast, region ±→genes). ✅
2. Does the contrast requirement extend to **TF ChIP / eCLIP** binding regulons? (pending)
3. **Batch 6** keep-or-drop once region-level lands (lean: drop, ship region-level).
4. **GitHub fork of all code** — final step, after the science is done.

## COMPLIANCE INVARIANTS
All shipped = NIH-funded + anonymously public; credential-walled (4DN/HuBMAP) & controlled-access excluded;
Tier-B (gnomAD/HPO) held; nothing derived from PIGEAN/EAGGL output (no circularity).
