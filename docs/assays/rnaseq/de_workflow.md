# RNA DE Workflow

`geneset-extractors workflows rna_de_prepare` is the upstream RNA differential-expression workflow.

It is not a new extractor. It prepares standardized differential-expression output that can be consumed directly by the existing extraction layer:

- `convert rna_deg`
- `convert rna_deg_multi`

The workflow writes:

- `deg_long.tsv`
- `comparison_manifest.tsv`
- `comparison_audit.tsv`
- `comparison_selected_samples.tsv`
- `prepare_summary.json`
- optional `pseudobulk_matrix.tsv`
- optional `pseudobulk_metadata.tsv`
- optional extractor outputs under `extractor/`

## Design boundary

The intended split is:

- `workflows/preprocessing`: start from counts plus metadata, fit DE, emit standardized long tables
- `convert`: start from already scored DE tables or precomputed program loadings, emit gene programs and GMT outputs

This keeps the extractor layer assay-result focused and keeps heavier upstream inference optional.

One important boundary now exists in two places:

- workflow-side `--de_mode` controls how DE is fit and which samples enter each contrast
- extractor-side `--postprocess_mode` on `rna_deg` / `rna_deg_multi` controls how an already-fit DE table is turned into signatures

They are intentionally separate. `--de_mode harmonizome` is not the same knob as `--postprocess_mode harmonizome`.

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

## DE modes

- `--de_mode modern` (default)
  - general-purpose behavior
  - uses all eligible samples for each contrast
  - allows covariates, batch columns, and repeated-measures designs subject to backend support
- `--de_mode harmonizome`
  - explicit bulk-RNA preset for GTEx/Harmonizome-style aging contrasts
  - balances each two-group comparison to `min(n_group_a, n_group_b)` after eligibility filtering
  - deterministic by seed; default seed is `1`
  - keeps the fit simple but still allows explicit fixed-effect covariates, for example `SEX,SMTSD`
  - warns if you use the preset without explicit covariates because broad tissues are more likely to drift toward generic signatures
  - writes selected sample IDs and pre/post balance counts to:
    - `comparison_selected_samples.tsv`
    - `comparison_audit.tsv`
  - currently limited to simple two-group bulk designs
  - rejects batch columns and repeated-measures flags so the fit stays auditable

Balancing is not the universal default. If you want the original general-purpose behavior, stay in `modern` mode or leave `--de_mode` unset.

Use `harmonizome` when the goal is conservative signature generation in broad, heterogeneous tissues where an imbalanced fit tends to surface mitochondrial, housekeeping, or generic bulk-expression signals. Do not use it as a universal default:

- it discards samples deliberately
- it is bulk-only in the current implementation
- it is not the right preset for repeated-measures or batch-heavy designs
- the general `modern` mode remains preferable when you want maximum power and standard DE inference

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

## Harmonizome-style bulk example

This mode exists to more closely match the published GTEx aging notebook behavior: broad-tissue stratified contrasts, balanced age-bin pools, and a simple two-group fit.

```bash
geneset-extractors workflows rna_de_prepare \
  --modality bulk \
  --counts_tsv path/to/tissue_counts.tsv \
  --sample_metadata_tsv path/to/sample_metadata.tsv \
  --sample_id_column sample_id \
  --feature_id_column gene_id \
  --group_column age_bin \
  --comparisons_tsv path/to/gtex_age_comparisons.tsv \
  --stratify_by tissue \
  --de_mode harmonizome \
  --covariates SEX,SMTSD \
  --balance_seed 1 \
  --backend auto \
  --out_dir results/rna_de_harmonizome \
  --organism human \
  --genome_build hg38
```

The audit outputs will record both the eligible pool and the balanced pool per contrast, along with the resolved covariates and the balancing seed.

A concrete validation example for the GTEx adipose aging contrast is recorded in:

- `docs/assays/rnaseq/harmonizome_validation_note.md`

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
  --covariates SEX,SMTSD \
  --backend auto \
  --out_dir results/rna_de_age_bins \
  --organism human \
  --genome_build hg38
```

This emits many contrast rows in one long table while preserving a generic interface. There is no GTEx-specific hardcoded preset in the extractor layer.

If you want GTEx-like contrasts but still want the general-purpose fit, keep `--de_mode modern` and use all eligible samples, optionally with explicit covariates such as `SEX,SMTSD`. If you want a narrower, more notebook-like and more conservative signature-oriented fit, switch to `--de_mode harmonizome`.

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
- use `--de_mode harmonizome` only when you explicitly want balanced notebook-style bulk contrasts
- if you use `harmonizome`, prefer explicit covariates for broad tissues rather than relying on the warning-only no-covariate path

scRNA:

- donor-level pseudobulk by default
- require meaningful donor counts per group
- avoid naive per-cell DE by default

Repeated measures:

- prefer `r_dream`
- otherwise fail clearly unless the user explicitly allows an approximate fallback
