# scRNA LIGER Preparation Workflow

This page documents the in-repo preparation path for deriving scRNA gene programs with LIGER/iNMF and then converting those loadings into DIG gene sets.

Command:

```bash
geneset-extractors workflows scrna_liger_prepare --help
```

## Purpose

`scrna_liger_prepare` standardizes the upstream preparation around LIGER:

- accepts one of `matrix_tsv`, `h5ad`, `seurat_rds`, or `mtx_dir`
- for `matrix_tsv`, performs deterministic downsampling and optional split-by-cell-type preparation in Python
- generates `run_liger.sh` scripts that call the packaged R/iNMF runner
- generates `run_geneset_extractors_from_liger.sh` scripts that feed `gene_loadings.tsv` into `convert rna_sc_programs`

The extractor-side `rna_sc_programs` converter remains ingestion-only.

## Required inputs

Choose exactly one primary expression input:

- `--matrix_tsv` plus `--meta_tsv`
- `--h5ad`
- `--seurat_rds`
- `--mtx_dir`

For `matrix_tsv`, recommended companion metadata columns are:

- `cell_id`
- `dataset` or `donor_id`
- `cell_type`

## Recommended matrix/metadata path

```bash
geneset-extractors workflows scrna_liger_prepare \
  --matrix_tsv path/to/cell_by_gene_logcounts.tsv \
  --meta_tsv path/to/cell_meta.tsv \
  --meta_cell_id_column cell_id \
  --dataset_column donor_id \
  --cell_type_column cell_type \
  --split_by_cell_type true \
  --max_cells_per_bucket 200 \
  --max_cells_total 50000 \
  --liger_k_grid 10,12,14,16,18,20,22,24 \
  --liger_top_n_genes 250 \
  --out_dir results/scrna_liger_prepare \
  --organism human \
  --genome_build hg38
```

Then run each generated subset script:

- `results/scrna_liger_prepare/subsets/<subset>/run_liger.sh`
- `results/scrna_liger_prepare/subsets/<subset>/run_geneset_extractors_from_liger.sh`

## Recommended direct h5ad path

```bash
geneset-extractors workflows scrna_liger_prepare \
  --h5ad path/to/cells.h5ad \
  --dataset_column donor_id \
  --cell_type_column cell_type__kp \
  --max_cells_total 50000 \
  --liger_top_n_genes 250 \
  --out_dir results/scrna_liger_prepare_h5ad \
  --organism human \
  --genome_build hg38
```

This mode creates a single `subsets/all/` runner and lets the R side split by cell type when the column is available.

## Output layout

- `prepare_summary.json`
- `subsets_manifest.tsv`
- `subsets/`
  - `cell_type=<CT>/counts_prefiltered.tsv` for `matrix_tsv` mode
  - `cell_type=<CT>/meta.tsv` for `matrix_tsv` mode
  - `run_liger.sh`
  - `run_geneset_extractors_from_liger.sh`
  - `liger_out/<cell_type>/gene_loadings.tsv` after execution
  - `liger_out/<cell_type>/gene_programs.txt` after execution
  - `liger_out/<cell_type>/cell_scores.tsv` after execution

## Conversion contract

The LIGER runner writes `gene_loadings.tsv` in a wide genes-by-program layout. That file is consumed via:

```bash
geneset-extractors convert rna_sc_programs \
  --liger_gene_loadings_tsv path/to/gene_loadings.tsv \
  --out_dir results/liger_gene_sets \
  --organism human \
  --genome_build hg38 \
  --score_transform positive \
  --select top_k \
  --top_k 250
```

## R dependencies

The generated `run_liger.sh` expects `Rscript` plus packages compatible with your current LIGER workflow, including:

- `rliger`
- `Seurat`
- `DelayedArray`
- `HDF5Array`
- `purrr`
- `dplyr`
- `clue`
- `proxy`
- `reticulate`
- `anndata`

## Notes

- `dataset_column` is the preferred generic batch field; `donor_id` can still be supplied there when that is the effective dataset identifier.
- `liger_top_n_genes` controls the size of the top-gene program exports written by the R runner.
- `extractor_top_k` controls the final DIG gene-set size emitted by `run_geneset_extractors_from_liger.sh`.
