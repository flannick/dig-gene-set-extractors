# HANDOFF_MANIFEST9 — Batch 9 (remaining clean gene sets; consolidates/retires Batch 3)

Date: 2026-07-01
Author: Gage (gage@broadinstitute.org)
Purpose: package the clean, not-yet-uploaded gene sets. **Retires Batch 3** — its GlyGen-observed + RNA-seq
parts are here; its histone parts (H3K9, histone×GTEx) are superseded by the region-level Batch 7.

## Contents
| Deliverable | Sets | Source (NIH, public) |
|---|---|---|
| GlyGen observed glycoproteins (+ N/O, ×GTEx) | 162 | GlyGen (NIGMS/Common Fund) + UniProt crosswalk + GTEx |
| ENCODE RNA-seq expressed genes (per biosample) | 296 | ENCODE RNA-seq (NHGRI) + NCBI gene_info |
| ENCODE RNA-seq consensus × GTEx (if present) | ~53 | ENCODE RNA-seq + GTEx |
| ENCODE-rE2G regulatory targets — background-corrected, symbol-cleaned | 508 | ENCODE-rE2G (NHGRI); GTEx symbol universe for cleanup |

## Notes
- **rE2G** here is the gene-level background-corrected version (Up/Down vs cross-experiment prevalence),
  with identifiers cleaned to valid HGNC symbols (original had ENSG/enhancer-element ID pollution). It is
  NOT region-level'd (rE2G is a regulatory-link model output, not raw accessibility peaks).
- **Retire/replace:** do not upload `batch3_derived_genesets.zip` — this batch supersedes it.

## Compliance
All ENCODE/NHGRI, GlyGen/NIGMS, GTEx/Common Fund, NCBI/NLM — anonymously public, NIH-funded. UniProt used
only as an ID crosswalk (interoperability). Nothing credential-walled or PIGEAN/EAGGL-derived.
