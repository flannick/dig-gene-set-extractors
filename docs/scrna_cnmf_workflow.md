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
  --cnmf_k_list auto \
  --cnmf_k auto \
  --seed 1 \
  --matrix_value_type logcounts \
  --out_dir results/scrna_cnmf_prepare
```

## K selection (default auto)

`K` is the number of cNMF components/programs. The recommended path is:

1. run `cnmf prepare/factorize/combine/k_selection_plot` over a K grid
2. select K from stability vs prediction-error tradeoff
3. run `cnmf consensus` at selected K

The default in this repo is reproducible auto-selection with:

- strategy: `largest_stable`
- threshold: `stability >= 0.95 * max(stability)`
- tie-break: largest K among candidates
- fallback: max stability when no candidate passes threshold

You can run selector directly:

```bash
omics2geneset workflows cnmf_select_k \
  --cnmf_output_dir results/scrna_cnmf_prepare/subsets/cell_type=<CT>/cnmf_out \
  --name <RUN_NAME>
```

Artifacts written in the cNMF run directory:

- `omics2geneset_k_selection.tsv`
- `omics2geneset_selected_k.json`

## Output layout

- `prepare_summary.json`
- `subsets_manifest.tsv`
- `subsets/`
  - `cell_type=<CT>/counts_prefiltered.tsv`
  - `cell_type=<CT>/meta.tsv`
  - `cell_type=<CT>/run_cnmf.sh`
  - `cell_type=<CT>/run_cnmf_consensus_auto_k.sh`
  - `cell_type=<CT>/run_omics2geneset_from_cnmf.sh`
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

Generated scripts:

- `run_cnmf.sh`: prepare/factorize/combine/k_selection_plot
- `run_cnmf_consensus_auto_k.sh`: selects K via `workflows cnmf_select_k` (or fixed `--cnmf_k`) and runs consensus
- `run_omics2geneset_from_cnmf.sh`: finds consensus spectra file (`tpm` by default) and runs `convert rna_sc_programs`

## Workflow flags and defaults

| Flag | Meaning | Default |
|---|---|---|
| `--cnmf_k_list` | K grid for prepare (`auto` uses subset-size rule) | `auto` |
| `--cnmf_k` | K for consensus script (`auto` or integer) | `auto` |
| `--cnmf_n_iter` | cNMF iterations for prepare | `100` |
| `--cnmf_numgenes` | Number of high-variance genes for cNMF | `2000` |
| `--seed` | Reproducible downsampling and cNMF seed passthrough | `1` |
| `--cnmf_total_workers` | Total workers for factorize command | `1` |
| `--cnmf_local_density_threshold` | Consensus local density threshold | `0.5` |
| `--cnmf_local_neighborhood_size` | Consensus local neighborhood size | `0.3` |

Auto K-grid by subset size (`n_cells_after_downsampling`):

- `<= 1500`: `5 8 10 12 15`
- `1500-5000`: `10 15 20 25 30`
- `> 5000`: `15 20 25 30 35 40 45 50`

Override examples:

```bash
# fixed grid + fixed K
omics2geneset workflows scrna_cnmf_prepare \
  ... \
  --cnmf_k_list "10 15 20 25" \
  --cnmf_k 20

# auto grid + stricter local-max K selection
omics2geneset workflows cnmf_select_k \
  --cnmf_output_dir <CNMF_OUTDIR> \
  --name <RUN_NAME> \
  --strategy largest_stable \
  --stability_frac_of_max 0.95 \
  --require_local_max true
```

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
