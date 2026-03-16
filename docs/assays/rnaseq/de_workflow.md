# RNA DE Workflow

`geneset-extractors workflows rna_de_prepare` is the upstream RNA differential-expression workflow.

It is not a new extractor. It prepares standardized differential-expression output that can be consumed directly by the existing extraction layer:

- `convert rna_deg`
- `convert rna_deg_multi`

The workflow writes:

- `deg_long.tsv`
- `comparison_manifest.tsv`
- `prepare_summary.json`
- optional `pseudobulk_matrix.tsv`
- optional `pseudobulk_metadata.tsv`
- optional extractor outputs under `extractor/`

## Design boundary

The intended split is:

- `workflows/preprocessing`: start from counts plus metadata, fit DE, emit standardized long tables
- `convert`: start from already scored DE tables or precomputed program loadings, emit gene programs and GMT outputs

This keeps the extractor layer assay-result focused and keeps heavier upstream inference optional.

## Supported v1 inputs

- `--modality bulk`
- `--modality scrna`
- dense TSV count matrices only
- TSV metadata only
- no AnnData/Scanpy requirement

## Standardized DE output

`deg_long.tsv` is written in long format with at least:

- `comparison_id`
- `gene_id`
- `gene_symbol`
- `logFC`
- `stat`
- `pvalue`
- `padj`

Helpful columns are also written when available:

- `group_a`
- `group_b`
- `stratum`
- `backend`
- `n_group_a`
- `n_group_b`
- `mean_expr`
- `model_formula`

## Backends

- `--backend auto`
  - prefers implemented/available R backends
  - otherwise falls back to `lightweight`
- `--backend lightweight`
  - approximate numpy-based fixed-effect backend
  - low-expression filtering
  - library-size normalization
  - log-scale modeling
  - BH correction within each emitted contrast
- `--backend r_limma_voom`
  - generates and executes a limma/voom script when the required R packages are available
- `--backend r_dream`
  - generates and executes a dream/variancePartition script when the required R packages are available

Repeated-measures designs do not silently fall back to the naive lightweight method. If only `lightweight` is available, the workflow fails unless:

```bash
--allow_approximate_repeated_measures true
```

## Lightweight backend limitations

The lightweight backend is useful for small, simple, dependency-light runs, but it is still approximate.

Limitations:

- fixed-effect only
- no true random-effects inference
- normal-approximation p-values
- intended for simple two-group or reference-level designs
- should not be treated as a replacement for limma-voom or dream on complex studies

Use `r_limma_voom` or `r_dream` when design complexity matters and the R stack is available.

## Bulk example: two-group DE

```bash
geneset-extractors workflows rna_de_prepare \
  --modality bulk \
  --counts_tsv path/to/bulk_counts.tsv \
  --sample_metadata_tsv path/to/sample_metadata.tsv \
  --sample_id_column sample_id \
  --feature_id_column gene_id \
  --group_column condition \
  --comparison_mode condition_a_vs_b \
  --condition_a case \
  --condition_b control \
  --backend auto \
  --out_dir results/rna_de_bulk \
  --organism human \
  --genome_build hg38
```

## scRNA example: donor-level pseudobulk within cell type

```bash
geneset-extractors workflows rna_de_prepare \
  --modality scrna \
  --counts_tsv path/to/cell_by_gene_counts.tsv \
  --cell_metadata_tsv path/to/cell_metadata.tsv \
  --cell_id_column cell_id \
  --feature_id_column gene_id \
  --donor_column donor_id \
  --cell_type_column cell_type \
  --group_column condition \
  --comparison_mode condition_a_vs_b \
  --condition_a treated \
  --condition_b control \
  --min_cells_per_pseudobulk 20 \
  --min_donors_per_group 3 \
  --backend auto \
  --out_dir results/rna_de_scrna \
  --organism human \
  --genome_build hg38
```

By default, this path aggregates counts to donor-level pseudobulks, optionally split by cell type when `--cell_type_column` is provided.

## GTEx-like example: age bins within tissue

```bash
geneset-extractors workflows rna_de_prepare \
  --modality bulk \
  --counts_tsv path/to/tissue_counts.tsv \
  --sample_metadata_tsv path/to/sample_metadata.tsv \
  --sample_id_column sample_id \
  --feature_id_column gene_id \
  --group_column age_bin \
  --comparison_mode reference_level \
  --reference_level 20-29 \
  --stratify_by tissue \
  --covariates sex \
  --backend auto \
  --out_dir results/rna_de_age_bins \
  --organism human \
  --genome_build hg38
```

This emits many contrast rows in one long table while preserving a generic interface. There is no GTEx-specific hardcoded preset in the extractor layer.

## Run the extractor after prepare

If you want the workflow to emit gene programs immediately, use:

```bash
geneset-extractors workflows rna_de_prepare \
  --modality bulk \
  --counts_tsv path/to/bulk_counts.tsv \
  --sample_metadata_tsv path/to/sample_metadata.tsv \
  --sample_id_column sample_id \
  --feature_id_column gene_id \
  --group_column condition \
  --comparison_mode condition_a_vs_b \
  --condition_a case \
  --condition_b control \
  --backend lightweight \
  --run_extractor true \
  --extractor_padj_max 0.05 \
  --extractor_min_abs_logfc 0.5 \
  --extractor_gmt_split_signed true \
  --extractor_gmt_topk_list 200 \
  --extractor_gmt_min_genes 100 \
  --extractor_gmt_max_genes 500 \
  --out_dir results/rna_de_bulk_with_extract \
  --organism human \
  --genome_build hg38
```

This writes the prepared DE artifacts first, then calls `convert rna_deg_multi` on the emitted `deg_long.tsv`.

## Extractor-side compatibility

The standardized long DE table is designed to work directly with `rna_deg_multi`.

The RNA DE converters now also support:

- `--padj_max`
- `--pvalue_max`
- `--min_abs_logfc`
- `--score_mode signed_neglog10padj`

Common column names are auto-detected more robustly, including:

- `logFC`
- `log2FoldChange`
- `avg_log2FC`
- `padj`
- `adj.P.Val`
- `FDR`

## Recommended defaults

Bulk:

- use all samples
- use explicit reference levels when needed
- keep low-expression filtering on
- prefer `auto` backend resolution

scRNA:

- donor-level pseudobulk by default
- require meaningful donor counts per group
- avoid naive per-cell DE by default

Repeated measures:

- prefer `r_dream`
- otherwise fail clearly unless the user explicitly allows an approximate fallback
