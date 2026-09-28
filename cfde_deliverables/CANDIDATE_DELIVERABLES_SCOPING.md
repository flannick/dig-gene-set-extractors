# Candidate Gene-Set Deliverables — Scoping Sheet

NIH-interoperability showcase. Rule: **NIH-funded + public only**; **don't recreate maintained published gene-set products** — but **deriving gene sets FROM public raw/processed data (cited) is allowed** when no usable gene-set format exists. Big data = **write pipeline here, run on cluster** (Ryan's pattern).

## Done / in hand
| Deliverable | Source | Funding | Status |
|---|---|---|---|
| Glyco **A/B/C/D** × 53 GTEx tissues | GTEx + UniProt + NCBI RefSeq + PIGEAN T2D | NIH (Common Fund / NLM / NIDDK) | ✅ built, zipped, dual-sourced |
| SenNet **SenMayo / SenSkin** standalone + × GTEx | SenNet/Saul/Wyles | NIH (Common Fund / NIA) | ✅ cited (published sets, not recreated) |
| **ENCODE ATAC accessible-genes** (K562/HepG2/HCT116) | ENCODE ATAC peaks + UCSC refGene hg38 | NHGRI/NIH | ✅ prototype + generalized extractor `accessible_genes/scripts/derive_accessible_genes.py`; scale = loop 463 biosamples on cluster |

## Candidate list (fed by user; all NIH + public; all cluster-bound)
| # | Source | Modality → extractor | Est. sets | Funding | Don't-redo check | Notes |
|---|---|---|---|---|---|---|
| — | **4DN-original** ATAC + HiCAR | accessibility → our `derive_accessible_genes.py` (+HiCAR anchors) | ~430 exps | NIH Common Fund | 4DN doesn't package as gene sets | verified 4DN-generated, not ENCODE-rehosted |
| 9 | **GEO** curated scRNA (20–30 Series) | programs → `scrna_cnmf_prepare`+`rna_sc_programs` (+`sc_rna_marker`) | 300–900 | GEO = NCBI/NLM/NIH infra, public; **per-study funding varies** | dedup vs atlas papers | heterogeneous supp-file formats = harder wrapper |
| 12 | **NeMO/BICCN** brain scRNA | programs → `scrna_cnmf_prepare`+`rna_sc_programs` | 500–1,500 | BRAIN Initiative (NIH) | BICCN publishes markers; cNMF programs = value-add | heaviest (cNMF compute) |
| 13 | **NeMO/BICCN** brain scATAC | accessibility → `atac_sc_10x`/`atac_bulk_matrix` (+ our peak→gene) | 200–800 | BRAIN Initiative (NIH) | check NeMO/BICCN | same modality as ENCODE work |
| 14 | **HuBMAP** scRNA/snRNA | programs → `scrna_cnmf_prepare`+`rna_sc_programs` | 500–2,000 | HuBMAP (NIH Common Fund) | HuBMAP publishes cell-type markers | 5,032 datasets, 27 organs |
| 15 | **HuBMAP** scATAC/ATAC | accessibility → `atac_sc_10x`/`atac_bulk`/`atac_bulk_matrix` (+ our peak→gene) | 200–800 | HuBMAP (NIH Common Fund) | check HuBMAP | GTF peak→gene = exactly our extractor |
| 16 | **HTAN** scRNA tumor-state | programs → `scrna_cnmf_prepare`+`rna_sc_programs` (+`sc_rna_marker`) | 200–700 | NCI/NIH (Cancer Moonshot) — **⚠️ TIERED ACCESS** | check HTAN cell-state annots | **use OPEN/Level-3+ tier ONLY; raw/Level-1-2 is dbGaP-CONTROLLED — do not ingest** (V7.0: 14 atlases, 20 sites, 2,372 cases, 31 assays) |

## Key efficiency insight
**One ATAC extractor covers all accessibility deliverables** — ENCODE, 4DN, NeMO #13, HuBMAP #15 differ only in peak source. `derive_accessible_genes.py` (promoter TSS±1kb peak overlap, GTF/refGene) is the shared core; scATAC variants just need pseudobulk + cell-type aggregation upstream (`atac_sc_10x`). So #13/#15 reuse most of what's built.

The **scRNA-program builds (#12, #14)** are the heavier, distinct effort (cNMF on single-cell matrices).

## Sequencing vs ~9-day showcase
- **Now/fast:** ENCODE + 4DN accessible-genes (extractor exists); HuBMAP/NeMO scATAC (#13/#15) reuse it once peaks pulled on cluster.
- **Cluster-run, larger:** scRNA programs (#12, #14) — write pipeline here, run on cluster, lands during/after showcase.
- **Every derived set's `provenance.json` must cite the source dataset (accession/UUID/DOI) + funder.**
