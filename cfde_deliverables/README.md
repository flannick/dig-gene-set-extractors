# CFDE gene-set extractors — derivation code

Code used to derive publicly-available, NIH-funded gene-set deliverables for the CFDE data-interoperability
showcase. Every deliverable is emitted in the **contract format** — one directory per set containing:

```
geneset.tsv            # header 'gene' + one HGNC symbol per line
genesets.gmt           # <name>\t<description>\t<gene1>\t<gene2>...
geneset.meta.json      # standard_name, library, description, n_genes, method, caveat, source, ...
geneset.provenance.json# operation, inputs, funding, source_citation, public:true
```

Sets are packaged into numbered `HANDOFF_MANIFEST*` batches (see `manifests/`) with md5 verification.

## Repository layout
| Path | What |
|---|---|
| `local/` | Non-cluster derivations (run on a laptop/workstation; mostly Python stdlib) |
| `cluster/` | Cluster pipelines: scRNA→cNMF gene programs, CATLAS full-222 scATAC accessibility |
| `reserve/` | Held/pending scripts (auth-gated or awaiting review) — **not shipped** |
| `manifests/` | `HANDOFF_MANIFEST*.md` — batch contents + compliance notes |
| `METHODS_accessibility_control.md` | The background-prevalence control method (shared design doc) |
| `PROJECT_STATUS.md` | Ledger: uploaded / held / running / dropped |
| `CANDIDATE_DELIVERABLES_SCOPING.md` | Candidate-deliverable list + in/out-of-scope rationale |

## Method highlights
- **Background-prevalence control** (accessibility, TF/eCLIP binding, CATLAS scATAC): a gene/region
  ubiquitously open/bound across experiments is treated as background. Per experiment, leave-one-out:
  `Up` = present here AND cross-experiment prevalence < 0.25 (specific); `Down` = absent here AND
  prevalence > 0.75. See `METHODS_accessibility_control.md`. **Honest scope: relative specificity, not a
  GC/library-normalized read-count differential.**
- **Region-level accessibility** kept as BEDs (peak ± regions) and mapped to genes via promoter overlap.
- **Interoperability crosses**: corrected regulons × disease/mito/dosage/glyco, chromatin states,
  eQTL/sQTL × accessibility convergence, master-regulator-per-tissue, multi-evidence disease genes.
- **scRNA gene programs**: self-contained cNMF (scanpy/anndata/cnmf) → per-program top genes.
- **CATLAS full-222**: fragment aggregation (GSE184462) → per-cell-type promoter accessibility →
  background control (see `cluster/cluster_catlas/README_catlas.md`).

## Script → deliverable map (selected)
| Script | Deliverable |
|---|---|
| `local/accessibility_bgcontrast/scripts/derive_accessibility_background_contrast.py` | gene-level ATAC/DNase/histone/rE2G Up/Down (corrected) |
| `local/accessibility_regions/scripts/derive_accessibility_regions.py` | region-level accessibility + BEDs (Batch 7) |
| `local/concordance_genesets/scripts/*` | accessibility × expression concordance (Batch 8) |
| `local/encode_regulons_bgcontrast/scripts/derive_regulon_background_contrast.py` | corrected TF/eCLIP regulons (Batch 10) |
| `local/interop_intersections/scripts/derive_*.py` | 10 interoperability crosses (Batch 11) |
| `local/catlas_genesets/scripts/derive_catlas_genesets.py` | CATLAS 28-subset (Batch 12) |
| `local/gtex_qtl`, `clinvar_genesets`, `clingen_dosage`, `glygen_genesets`, `mito_genesets` | GTEx/ClinVar/ClinGen/GlyGen/MitoCarta sets (Batch 4/9) |
| `cluster/cluster_pipeline/scrna_cnmf_programs.py` | BICCN scRNA cNMF programs |
| `cluster/cluster_catlas/*` | CATLAS full-222 scATAC accessibility |

## Reproducing
Local scripts fetch their own public inputs (ENCODE/GTEx/ClinVar/UCSC/NCBI/GlyGen) and write a
`<deliverable>/output/` tree, then are zipped into a handoff batch. See each manifest for exact contents.
Cluster pipelines: see the `README`/`setup_env.sh` under `cluster/`.

## Compliance
All inputs are **NIH-funded and anonymously public** (no credential-walled data). Deliverables are newly
derived here (not re-releases of published gene-set products) and are **not** PIGEAN/EAGGL-derived.
Peak/binding-based sets use the background control above. See per-manifest compliance notes.
