# RESERVE — hold for later (NOT for the incubator / NOT for upload today)

Date: 2026-06-30
Decision (gage): these are **kept for later**, deliberately excluded from the NIH incubator/showcase
handoffs and from today's uploads. Do **not** package, zip, or upload anything in this directory without
explicit re-authorization.

## Contents

### `iptmnet_genesets/`
PTM enzyme/effector gene sets derived from iPTMnet (15 PTM-type substrate sets, per-enzyme sets,
enzyme/effector standalone, ×GTEx). **Held pending license/provenance sign-off** — iPTMnet is
CC BY-NC-SA 4.0 with HPRD-dominant source content (restrictive). Will not be uploaded today or for the
incubator. See `iptmnet_genesets/LICENSE_REVIEW.md`, `TEAM_BLURB.md`.

### `encode_crispr_kd_regulons/`
ENCODE **CRISPR + CRISPRi** RNA-seq knockdown regulons (K562/HepG2 cell-line perturbation; direct+indirect
effects). Same engine as the shippable shRNA regulons, but the non-shRNA assays are **held for later** per
gage. Public ENCODE/NHGRI data; the hold is a scoping decision, not a licensing one. Fully labeled
`evidence_type: experimental_perturbation`, `biosample: K562/HepG2`.

### `gencode_antisense/`
GENCODE cis-antisense → protein-coding target map (ncRNA regulatory relationships) + ×GTEx. Public
GENCODE/NHGRI data. **Held with the iPTM/non-shRNA reserve** as "ncRNA stuff for later" per gage.

## What ships instead (in the main project, incubator-bound)
shRNA RNA-seq KD regulons (labeled), ENCODE TF-ChIP regulons, eCLIP RBP regulons, histone×GTEx,
RNA-seq×GTEx, GlyGen observed glycoproteins, MitoCarta-derived mito sets — all NIH-funded, public,
tissue/organism-grounded or explicitly labeled.

---

## Added 2026-06-30 (compliance holds)

### `tier_b_pending_review/`
Sets derived from **NIH-CO-funded / multi-funder open resources** — held pending gage's funding-scope
sign-off (likely fine for CFDE, holding to be safe). Contains:
- `gnomad_constraint/` (gnomAD LoF/missense constraint sets)
- `hpo_genesets/` (5,261 HPO phenotype gene sets)
- `interop_intersections/` (the gnomAD-dependent intersections, incl. constraint×disease/QTL/tissue + the disease∩constrained∩tissue triples, and GlyGen×gnomAD)
NOTE: **UniProt is NOT held** — it is used only as an ID *crosswalk* (AC→gene), an interoperability tool, so GlyGen-observed and other crosswalk-using sets stay shippable.

### `cluster_auth_gated/`
- `run_4dn_atac.sh` — 4DN ATAC pipeline. 4DN downloads are **access-key gated (403 anon, confirmed)**; gage lacks authorization → held until authorized or rewritten to the released-only open S3.
- **HuBMAP (reserved, no local files):** HuBMAP processed assets are **403 anon** (Globus/auth) → not anonymously public; excluded from the scRNA wiring, revisit shortly.

### scRNA datasets that ARE clean (wired, not reserved)
Allen/BICCN human **M1 10x** and human **cortex SMART-seq** — anon-public on the Allen S3 (verified),
NIH BRAIN/BICCN-funded. Wired as input-specified workflows in `cluster_pipeline/`.
