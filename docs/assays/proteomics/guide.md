# Proteomics Guide

This assay family now has two distinct entrypoints:

- `proteomics_diff`: legacy gene-level protein abundance table -> gene weights
- `ptm_site_diff`: site-level PTM contrast table -> compact gene program and GMT sets
- `ptm_site_matrix`: site-by-sample PTM matrix + sample metadata -> one or more contrast-specific gene programs

`ptm_site_diff` is the recommended PTM entrypoint when you already have a site-level contrast table. It is designed for phosphoproteomics first, but it also supports other site-resolved PTM tables when the input includes stable site identifiers or enough fields to reconstruct them.

`ptm_site_matrix` is the recommended PTM entrypoint for public site-by-sample matrices. It estimates within-study site-level contrasts, then reuses the same downstream scorer as `ptm_site_diff`.

`workflows ptm_prepare_public` is the workflow-layer normalizer for official CDAP/PDC public phosphosite and proteome reports. It converts those public reports into the standardized `ptm_matrix.tsv` / `sample_metadata.tsv` / optional `protein_matrix.tsv` interface consumed by `ptm_site_matrix`.

## Quickstart

### 1. Optional: build a local phosphosite bundle

The PTM bundle is optional. It improves harmonization and site ubiquity weighting when local prior tables are available.

```bash
geneset-extractors workflows ptm_prepare_reference_bundle \
  --sources_tsv <sources.tsv> \
  --out_dir <ptm_bundle_dir> \
  --organism human \
  --ptm_type phospho \
  --bundle_id phosphoproteomics_human_v1
```

This writes compact files such as:

- `phosphosite_aliases_human_v1.tsv.gz`
- `phosphosite_ubiquity_human_v1.tsv.gz`
- `bundle_provenance.json`
- `local_resources_manifest.json`

If you point `--resources_dir` at that directory, both `ptm_site_diff` and `ptm_site_matrix` will auto-resolve the two PTM priors by default for human phosphoproteomics.

### 2. Convert a site-level PTM differential table

```bash
geneset-extractors convert ptm_site_diff \
  --ptm_tsv <ptm_site_diff.tsv> \
  --out_dir <out_dir> \
  --organism human \
  --genome_build human \
  --ptm_type phospho \
  --resources_dir <ptm_bundle_dir>
```

Minimal no-bundle run:

```bash
geneset-extractors convert ptm_site_diff \
  --ptm_tsv <ptm_site_diff.tsv> \
  --out_dir <out_dir> \
  --organism human \
  --genome_build human \
  --ptm_type phospho \
  --use_reference_bundle false
```

### 2b. Standardize public CDAP/PDC PTM reports

Explicit local files:

```bash
geneset-extractors workflows ptm_prepare_public \
  --input_mode cdap_files \
  --ptm_report_tsv <Phosphoproteome.phosphosite.tmt11.tsv> \
  --protein_report_tsv <Proteome.tmt11.tsv> \
  --sample_design_tsv <study.sample.txt> \
  --sample_annotations_tsv <sample_annotations.tsv> \
  --out_dir <prepared_dir> \
  --organism human \
  --ptm_type phospho \
  --study_id <study_id> \
  --study_label <study_label>
```

Local PDC-manifest discovery mode:

```bash
geneset-extractors workflows ptm_prepare_public \
  --input_mode pdc_manifest \
  --pdc_manifest_tsv <manifest.tsv> \
  --source_dir <download_dir> \
  --sample_annotations_tsv <sample_annotations.tsv> \
  --out_dir <prepared_dir> \
  --organism human \
  --ptm_type phospho \
  --study_id <study_id>
```

This workflow writes:

- `ptm_matrix.tsv`
- `sample_metadata.tsv`
- optional `protein_matrix.tsv`
- `sample_id_map.tsv`
- `site_id_map.tsv`
- `prepare_summary.json`
- `bundle_source_row.tsv`

`bundle_source_row.tsv` is deliberately shaped so you can concatenate many prepared studies and feed them into `workflows ptm_prepare_reference_bundle`.

This workflow also performs assay-type QC for phosphoproteomics staging. `prepare_summary.json` reports:

- `dominant_residue_family`
- `fraction_s`
- `fraction_t`
- `fraction_y`
- `fraction_k`
- `phospho_like_fraction`

Default behavior is `--assay_type_policy warn`. For stricter public-data screening, use:

```bash
geneset-extractors workflows ptm_prepare_public \
  ... \
  --assay_type_policy fail
```

This is mainly intended to catch lysine-dominant or otherwise non-phospho-like public reports before they enter the phospho extractor path.

### 2c. Convert a standardized site-by-sample PTM matrix

Single study-wide condition contrast:

```bash
geneset-extractors convert ptm_site_matrix \
  --ptm_matrix_tsv <prepared_dir>/ptm_matrix.tsv \
  --sample_metadata_tsv <prepared_dir>/sample_metadata.tsv \
  --protein_matrix_tsv <prepared_dir>/protein_matrix.tsv \
  --study_contrast condition_a_vs_b \
  --condition_a case \
  --condition_b control \
  --out_dir <out_dir> \
  --organism human \
  --genome_build human \
  --ptm_type phospho \
  --resources_dir <ptm_bundle_dir>
```

When `--protein_matrix_tsv` is present, the default matrix behavior is:

- `--protein_adjustment_run_mode compare_if_protein`

This emits paired variants for the same biological contrast:

- `protein_adjustment=none`
- `protein_adjustment=subtract`

If you want a single adjusted or unadjusted run product instead, set:

```bash
--protein_adjustment_run_mode single
```

Per-group condition contrasts:

```bash
geneset-extractors convert ptm_site_matrix \
  --ptm_matrix_tsv <ptm_matrix.tsv> \
  --sample_metadata_tsv <sample_metadata.tsv> \
  --study_contrast condition_within_group \
  --group_column group \
  --condition_column condition \
  --condition_a case \
  --condition_b control \
  --out_dir <out_dir> \
  --organism human \
  --genome_build human \
  --ptm_type phospho
```

### 3. Validate outputs

```bash
geneset-extractors validate <out_dir>
```

## Input contract for `ptm_site_diff`

The extractor expects a TSV or TSV.GZ with site-level PTM rows and enough information to do three things:

1. construct or resolve a canonical site key
2. compute a site score
3. assign the site to one or more genes

Common columns include:

- `site_id`
- `site_group_id`
- `gene_id`
- `gene_symbol`
- `protein_accession`
- `residue`
- `position`
- `log2fc`
- `stat`
- `pvalue` or `padj`
- `localization_prob`
- `peptide_count`
- `protein_log2fc` or `protein_stat`

The extractor will use explicit column flags when provided, otherwise it falls back to common defaults.

## Input contract for `ptm_site_matrix`

`ptm_site_matrix` expects:

1. a wide PTM matrix TSV with one row per site and one column per sample
2. a sample metadata TSV with at least `sample_id`, and usually `group` and/or `condition`
3. optionally, a matched protein matrix with one row per protein and the same sample columns

The PTM matrix should contain the same site metadata used by `ptm_site_diff`, for example:

- `site_id`
- `gene_id`
- `gene_symbol`
- `protein_accession`
- `residue`
- `position`
- optional confidence fields such as `localization_prob` and `peptide_count`

The sample metadata should contain:

- `sample_id`
- `condition` for `condition_a_vs_b` or `condition_within_group`
- `group` for `group_vs_rest` or `condition_within_group`

The current matrix front-end supports:

- `--study_contrast condition_a_vs_b`
- `--study_contrast group_vs_rest`
- `--study_contrast condition_within_group`
- `--study_contrast baseline`

For matrix contrasts, the default effect metric is `welch_t`; `mean_diff` is available when you want a simpler within-study contrast.

## Public-report workflow contract

`ptm_prepare_public` is intentionally narrow. It is for official local CDAP/PDC public reports, not for already standardized matrices.

Inputs supported in v1:

- phosphosite report TSV
- optional matched proteome report TSV
- optional `*.sample.txt` experimental design table
- optional sample-annotation sidecar
- optional PDC manifest plus `--source_dir` for local discovery

The helper does structural harmonization only:

- parse `Phosphosite` strings into `site_id` or `site_group_id` when lossless
- preserve `raw_site_label`, `raw_peptide`, and `site_parser_status`
- normalize the public reports into repo-standard wide matrices
- enrich sample metadata from `sample.txt` and optional sidecars
- keep the source numeric values unchanged

The helper does not do:

- cross-study alias resolution
- site ubiquity weighting
- protein-adjusted PTM scoring
- duplicate collapse
- gene aggregation

Those remain the responsibility of `ptm_site_diff`, `ptm_site_matrix`, and the optional PTM bundle.

## Recommended defaults

The default PTM configuration is intended to produce connectable, direction-aware phosphoregulation rather than a raw count of changed sites:

- `--score_mode auto`
- `--score_transform signed`
- `--protein_adjustment subtract`
- `--protein_adjustment_run_mode compare_if_protein`
- `--protein_adjustment_lambda 1.0`
- `--confidence_weight_mode combined`
- `--site_dup_policy highest_confidence`
- `--gene_aggregation signed_topk_mean`
- `--gene_topk_sites 3`
- `--select top_k --top_k 200`
- `--normalize within_set_l1`
- `--emit_full true`
- `--emit_gmt true`
- `--gmt_split_signed true`
- `--use_reference_bundle true`

If the PTM bundle is missing under `--resource_policy skip`, the converter warns and falls back to direct canonicalization and no ubiquity penalty.

For site-count benchmarking without changing the main default, `ptm_site_matrix` also supports:

- `--emit_gene_topk_site_comparison true`
- `--gene_topk_site_compare_to 1`

This emits an additional stricter `gene_topk_sites=1` variant beside the default `gene_topk_sites=3`.

## Output contract

`ptm_site_diff` writes:

- `geneset.tsv`
- `geneset.meta.json`

Optional or default extra outputs:

- `geneset.full.tsv`
- `genesets.gmt`
- `run_summary.json`
- `run_summary.txt`

`geneset.tsv` contains the selected compact program with at least:

- `gene_id`
- `gene_symbol`
- `score`
- `weight`
- `rank`
- `n_supporting_sites`
- `top_sites`

`ptm_site_matrix` writes the same per-contrast files. When the study design produces multiple contrasts, the root output also contains:

- `manifest.tsv`
- `contrast_qc.tsv`
- root `genesets.gmt` aggregating emitted child GMTs
- root `run_summary.json` / `run_summary.txt`

When paired protein-adjustment or stricter site-cap variants are emitted, `manifest.tsv` and `contrast_qc.tsv` also include:

- `protein_adjustment`
- `gene_topk_sites`
- `variant_id`

## Conceptual model

A site-level PTM table is not a gene program by itself. The extractor builds one in the following stages:

1. compute a base site statistic
2. optionally weight it by confidence and localization
3. optionally subtract matched protein change
4. optionally downweight ubiquitous sites using a local phosphosite prior
5. collapse duplicate rows to canonical sites
6. aggregate canonical sites to genes with a site-count-robust rule
7. optionally compare unadjusted and protein-adjusted variants
8. select and normalize a compact gene program
9. optionally emit signed GMT sets

The most important design choice is the default site-to-gene aggregation. `signed_topk_mean` avoids overweighting proteins with many measured sites while still preserving coherent multi-site regulation.

Theory and equations: `docs/assays/proteomics/methods.tex`

## Resource bundle notes

The optional PTM bundle does two things in v1:

- harmonize heterogeneous site identifiers
- supply a detection-frequency style ubiquity prior

The recommended human phosphoproteomics resource ids are:

- `phosphosite_aliases_human_v1`
- `phosphosite_ubiquity_human_v1`

Preset:

- `phosphoproteomics_default_optional_human`

Bundle guide: `docs/assays/proteomics/reference_bundle.md`

## Common pitfalls

- Ambiguous gene assignments are dropped by default.
  - If your source table contains many multi-gene site groups, inspect `n_rows_ambiguous_gene` and consider `--ambiguous_gene_policy split_equal` only if that is biologically justified.
- Protein adjustment changes interpretation.
  - `--protein_adjustment subtract` makes the default output closer to PTM-specific regulation beyond total protein abundance change.
- Matrix mode defaults to paired adjusted vs unadjusted output when a matched protein matrix is available.
  - Set `--protein_adjustment_run_mode single` if you only want one variant.
- Missing bundle resources are non-fatal by default.
  - Under `--resource_policy skip`, the run continues and records missing resources in metadata/run summaries.
- GMT emission still has size guardrails.
  - Small signed subsets may be skipped unless you set `--emit_small_gene_sets true` or lower `--gmt_min_genes`.
- `proteomics_diff` is not a PTM substitute.
  - It is a legacy gene-level abundance converter and does not implement site harmonization, protein adjustment, or site ubiquity weighting.
- PTM outputs can still be dominated by sample composition rather than tumor-intrinsic phosphoregulation.
  - `run_summary.json` and `run_summary.txt` now report `composition_qc`, `composition_warning`, and `tumor_intrinsic_confidence`.
  - These are advisory QC flags, not hard filters.

## Common PTM matrix pitfalls

- Matrix columns must match `sample_id` values in the metadata.
  - Sample metadata rows that do not appear as PTM matrix columns are ignored.
- `condition_within_group` needs both `group` and `condition`.
  - By default the converter looks for columns literally named `group` and `condition`.
- Missing values can suppress site contrasts before PTM scoring even starts.
  - `--missing_value_policy min_present` is the easier default.
  - `--missing_value_policy drop` is stricter and will skip any site with missing values in compared samples.
- Matched protein adjustment only works if the protein matrix can be keyed back to the PTM sites.
  - Include `protein_accession`, `gene_id`, or `gene_symbol` in the protein matrix rows.
- Public PTM reports are not all phosphoproteomics assays.
  - `ptm_prepare_public` can warn or fail early on lysine-dominant or otherwise non-phospho-like studies via `--assay_type_policy`.

## Standard public-data path

The intended public phosphoproteomics path is now:

1. `geneset-extractors workflows ptm_prepare_public`
2. optional `geneset-extractors workflows ptm_prepare_reference_bundle`
3. `geneset-extractors convert ptm_site_matrix`

Use `ptm_site_diff` only when you already have a site-level contrast table rather than a public site-by-sample matrix.
