# HANDOFF_MANIFEST10 — Batch 10 (background-corrected TF ChIP + eCLIP regulons)

Date: 2026-07-01
Author: Gage (gage@broadinstitute.org)
Status: applies the SAME background-prevalence correction as ATAC/DNase to the binding regulons.
Replaces the held raw TF/eCLIP regulons (which were over-inclusive due to promiscuous/HOT-region binding).

## Method
Background = per-gene **binding prevalence across all factors** of the assay (how many TFs/RBPs target the
gene); a gene bound by most factors is a promiscuous/HOT-region artifact (binding analog of a ubiquitously-
open promoter). Leave-one-out. Per factor:
- **specific_Up** = bound by this factor AND prevalence(excl self) < 0.25 → specifically-bound targets
- **promiscuous_Down** = NOT bound here AND prevalence > 0.75 → commonly-bound-elsewhere genes

## Contents — 2,440 gene sets
| Assay | specific_Up | promiscuous_Down | promiscuous(>0.75) genes |
|---|---|---|---|
| TF ChIP-seq (1,149 factors) | 1,146 | 958 | 493 |
| eCLIP (168 RBPs) | 168 | 168 | 6,372 |

## Notes
- Same correction/rationale as the accessibility work (see METHODS_accessibility_control.md); relative
  binding *specificity*, not a read-count differential.
- The raw (uncorrected) TF/eCLIP regulons + their ×GTEx are retained locally but NOT shipped (superseded here).

## Compliance
ENCODE TF ChIP-seq + eCLIP (NHGRI), anonymously public, NIH-funded. No credential-walled data.
