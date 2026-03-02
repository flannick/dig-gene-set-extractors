# scRNA cNMF Preparation Workflow

This page documents the first-class in-repo preprocessing path for large scRNA atlases before cNMF.

Command:

```bash
omics2geneset workflows scrna_cnmf_prepare --help
```

## Purpose

`scrna_cnmf_prepare` prepares cNMF-ready subset matrices from a dense cell x gene TSV plus metadata:

- deterministic downsampling (seeded)
- optional split-by-cell-type runs
- donor-aware bucket balancing
- removal of zero-total cells and zero-total genes (required by cNMF best practice)
- per-subset `run_cnmf.sh` generation

`rna_sc_programs` remains ingestion-only and consumes the resulting cNMF gene spectra files.

## Required inputs

- `--matrix_tsv`: dense cell x gene matrix TSV.
  - header required
  - one row per cell
  - one gene column per feature
- `--meta_tsv`: metadata TSV with at least cell IDs.

## Recommended command

```bash
omics2geneset workflows scrna_cnmf_prepare \
  --matrix_tsv path/to/cell_by_gene_logcounts.tsv \
  --meta_tsv path/to/cell_meta.tsv \
  --meta_cell_id_column cell_id \
  --cell_type_column cell_type \
  --donor_column donor_id \
  --split_by_cell_type true \
  --max_cells_per_bucket 200 \
  --max_cells_total 20000 \
  --seed 1 \
  --matrix_value_type logcounts \
  --out_dir results/scrna_cnmf_prepare
```

## Output layout

- `prepare_summary.json`
- `subsets_manifest.tsv`
- `subsets/`
  - `cell_type=<CT>/counts_prefiltered.tsv`
  - `cell_type=<CT>/meta.tsv`
  - `cell_type=<CT>/run_cnmf.sh`
  - or `all/` when not splitting by cell type

## cNMF script generation

Each subset gets `run_cnmf.sh` with:

1. `cnmf prepare`
2. `cnmf factorize`
3. `cnmf combine`
4. `cnmf k_selection_plot`
5. commented `cnmf consensus` template
6. commented `omics2geneset convert rna_sc_programs ...` template

This repo does not require `cnmf` as a dependency. By default scripts are generated only.

Optional:

- `--execute true` attempts to execute generated scripts and fails early if `cnmf` is missing.

## Key tuning parameters

- `--split_by_cell_type`
- `--min_cells_per_cell_type`
- `--max_cells_per_bucket`
- `--max_cells_total`
- `--bucket_columns`
- `--seed`
- `--min_total_per_cell`
- `--min_total_per_gene`

## Notes on counts vs logcounts

- cNMF is naturally defined on nonnegative expression matrices.
- Counts are preferred when available.
- Logcounts can still be used for practical workflows if nonnegative.
- `scrna_cnmf_prepare` uses `--matrix_value_type` only for threshold defaults and warnings; it does not invert log transforms.
