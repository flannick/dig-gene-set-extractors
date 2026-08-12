# CATLAS full-222 cell-type accessibility gene sets (cluster)

Fragment-aggregation pipeline that builds background-corrected, per-cell-type promoter-accessibility
gene sets for **all 222 cell types** of the human single-cell chromatin-accessibility atlas
(Zhang et al. 2021 *Cell*; **GSE184462 / CATlas**; NIH-funded, anonymously public).

This supersedes the 28-cell-type subset (Batch 12) — same control method, full atlas.

## Data (downloaded on the cluster)
- `GSE184462_RAW.tar` (35.6 GB) — 155 per-sample fragment BEDs (`chrom start end barcode count`).
- `GSE184462_metadata.tsv.gz` (34 MB) — barcode→cell-type crosswalk
  (`cellID=sample+barcode | sample | replicate | logUMI | tsse | tissue | cell type | Life stage`;
  1,323,041 cells, 222 cell types, 155 samples).

## Join (verified)
- Sample: fragment file `GSM#####_<sample>_repN_fragments.bed.gz` → metadata sample `<sample>_N`
  (strip `GSM#####_`, `repN`→`_N`). Confirmed e.g. `adipose_omentum_SM-ADYHB_rep1` → `..._1` (3,799 cells).
- Barcode: fragment col4 == the part of `cellID` after `+` (same 22 bp). Confirmed matching prefixes.

## Method (same background control as ENCODE ATAC/DNase)
1. `catlas_make_promoters.py` — hg38 refGene → promoter windows (TSS ± `PROMOTER_WIN`, default 1 kb) →
   1 kb bin → gene map (`chrom:bin \t gene1;gene2`). Primary chromosomes only.
2. `catlas_aggregate_sample.sh` (per sample; UGER array) — for each fragment: barcode→cell type,
   position→promoter bin→gene, count **1 per fragment** overlapping a promoter. Emits
   `cell type \t gene \t count`.
3. `catlas_combine_call.py` — sum across samples; per-cell signal `= count / n_cells(cell type)`;
   **accessible** if `count ≥ MINCOUNT` (default 10) and `signal ≥ MINSIG` (default 0.15 frags/cell).
   Background = **cross-cell-type accessibility prevalence** (leave-one-out):
   - `CATLAS_<ct>_accessible_raw`  — accessible genes, **no control** (labeled raw)
   - `CATLAS_<ct>_specific_Up`     — accessible AND prevalence(excl self) < `LOW` (0.25): cell-type-specific
   - `CATLAS_<ct>_closed_Down`     — broadly-open (prev > `HIGH` 0.75) but closed here
   Cell types with `< MIN_CELLS` (25) cells dropped.

**Honest scope:** relative cell-type accessibility *specificity* (promoter fragment CPM/cell); **not** a
GC/library-normalized read-count differential. Same caveat as all our accessibility deliverables.

## Run
```bash
cd /humgen/diabetes2/users/gage/CFDE/cluster_catlas   # wherever these live
bash run_catlas_full.sh setup      # download 35.6 GB + metadata, untar, build promoters, filelist
bash run_catlas_full.sh run        # UGER array 1-155 (blocks via qsub -sync y); else local NPROC
bash run_catlas_full.sh combine    # -> gene sets -> zip -> copy to SDG + md5sum
# or: bash run_catlas_full.sh all
```
Tunable via env: `WORK`, `SDG`, `NPROC`, `MINSIG`, `MINCOUNT`, `LOW`, `HIGH`, `MIN_CELLS`, `PROMOTER_WIN`.

Expected output: `batch_catlas_full222.zip` = raw + specific_Up + closed_Down per cell type
(~3 × ~200 usable cell types), copied to the SDG submissions folder with an md5.

## Compliance
CATlas / GSE184462 — NIH-funded, anonymously public GEO. No credential-walled data, not PIGEAN/EAGGL-derived.
Validated locally on a synthetic mini-dataset (join, binning, multi-gene bins, per-cell normalization,
prevalence control, contract files) before the full run.
