# METHODS — Controlling chromatin accessibility (ENCODE ATAC/DNase) for gene-set derivation

_Project: CFDE gene-set deliverables. Audience: another Claude CLI / collaborator reviewing or reproducing
our accessibility handling. Status: advisor-approved design; region-level in production 2026-07-01._

## 1. The problem
ENCODE ATAC-seq and DNase-seq have **no matched input/genomic-DNA control** (unlike ChIP-seq). A naïve
gene set = "gene whose promoter overlaps a called peak" is therefore **over-inclusive**: promoters of
broadly-expressed / housekeeping genes are open in nearly every cell type, so the set is dominated by
ubiquitous accessibility and carries little cell-type-specific signal. Peaks alone are not a control.

## 2. What we do — an aggregated cross-experiment background contrast
We define accessibility **relative to an aggregated background built from all ENCODE experiments of the same
assay**, and emit matched **Up / Down** gene sets per experiment:

- **Background** = per-region (or per-gene-promoter) **accessibility prevalence** = fraction of ENCODE
  experiments (of that assay/mark) in which it is called accessible.
- **Leave-one-out**: a sample is excluded from its own background (avoids self-referential circularity).
- **Up** = accessible in this experiment AND prevalence < LOW (default 0.25) → *specifically accessible*.
- **Down** = NOT accessible here AND prevalence > HIGH (default 0.75) → *specifically closed where usually open*.
- **Histone ChIP** marks use a **per-mark** background (H3K4me3 vs H3K4me3 only), never pooled across marks.
- Thresholds are recorded in every set's `meta.json`.

This removes constitutively-open promoters (high prevalence → excluded from Up) — the de-biasing goal.

## 3. Two resolutions
1. **Gene-promoter level** (`accessibility_bgcontrast/`): background prevalence computed over genes
   (promoter TSS±1kb). Fast; derived from the per-biosample accessible-gene lists. (Batch 6.)
2. **Region level** (`accessibility_regions/`): re-download ENCODE peak BEDs (**kept** as reusable
   artifacts), bin the genome (1 kb), call each **region ± vs background**, then **map regions → genes**
   by promoter overlap. Writes per-experiment ±region BEDs + Up/Down gene sets. This is the rigorous,
   advisor-endorsed "region ± → map genes" design; the kept region BEDs support reuse in other projects.

## 4. HONEST SCOPE / limitations (important)
- This is a **peak-call prevalence background** = *relative accessibility specificity*. It is **NOT** a
  GC-content-aware, library-size-normalized, replicate-based **read-count differential** (DESeq2/edgeR),
  which the literature treats as the gold standard but which requires **raw reads/fragments** we do not
  re-process here. Every set is labeled with this caveat.
- The literature-sanctioned framings we align with: chromVAR's GC/accessibility-matched *background peak
  set* (Schep 2017), GREAT's region-*universe* for enrichment (McLean 2010), and MACS2's within-sample
  local-λ background (Zhang 2008). A single pooled track used as an *input-DNA proxy* is NOT established and
  risks circularity + washing out condition-specific peaks — we avoid that (leave-one-out; relative, not proxy).
- Region mapping is promoter-proximal; distal enhancer→gene (GREAT/Cicero/ABC) is a future upgrade.

## 5. Concordance extension (accessibility × expression)
Per matched ENCODE biosample we also compare accessibility to expression (ENCODE RNA-seq, same biosample):
`concordant_active` (open+expressed), `open_but_silent` (poised), `expressed_but_closed` (distal/latency).
Run in **two modes** — **raw** (uncontrolled accessibility) and **control** (background-corrected) — and a
**raw−control "background-sneak" delta** that catches accessibility signal that is actually background
(ubiquitous open) vs real signal absorbed by the control. Dir: `concordance_genesets/`.

## 6. Parameters
`LOW=0.25`, `HIGH=0.75` (prevalence thresholds), `BIN=1000` bp, promoter window `TSS±1000` bp, leave-one-out
on. All overridable via env; all recorded in provenance.

## 7. Code (all pure-Python stdlib, portable)
- `accessibility_bgcontrast/scripts/derive_accessibility_background_contrast.py` — gene-level contrast (GROUP_KEY for per-mark; SYM_UNIVERSE to clean rE2G identifiers)
- `accessibility_regions/scripts/derive_accessibility_regions.py` — region-level; keeps peak BEDs + ±region BEDs; atomic downloads; precomputed high-prevalence bins
- `accessibility_regions/scripts/run_histone_region.sh` — per-mark histone driver
- `concordance_genesets/scripts/derive_accessibility_expression_concordance.py` (raw+control)
- `concordance_genesets/scripts/derive_concordance_background_delta.py` (background-sneak)

## 8. Compliance
All inputs ENCODE/NHGRI (+ GTEx/Common Fund, NCBI/NLM), anonymously public, NIH-funded. No credential-walled
(4DN/HuBMAP) or controlled-access data. Nothing derived from PIGEAN/EAGGL output (these sets are inputs to
that pipeline — no circularity).

## 9. Reproduce (example)
```
# gene-level, ATAC:
IN_ZIP=ENCODE_ATAC_accessible_genes_20260630.zip LIB=ENCODE_ATAC_accessible_bgcontrast ASSAY="ATAC-seq" \
  python accessibility_bgcontrast/scripts/derive_accessibility_background_contrast.py
# region-level, ATAC (keeps BEDs):
ENC_ASSAY="ATAC-seq" ENC_OUTPUT="IDR thresholded peaks" LIB=ENCODE_ATAC_region \
  REFGENE=/path/refGene_hg38.txt.gz python accessibility_regions/scripts/derive_accessibility_regions.py
```
