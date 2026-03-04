# Drug Response Screens to Gene Sets (`drug_response_screen`)

This guide covers conversion of drug response screens (PRISM, GDSC/CTRP-like long tables)
into gene-weighted programs and GMT gene sets.

`drug_response_screen` builds gene programs from:

- response values per sample and drug
- drug-target annotations

and emits grouped outputs:

- `manifest.tsv`
- `program=<PROGRAM>/geneset.tsv`
- `program=<PROGRAM>/geneset.meta.json`
- optional per-program `geneset.full.tsv`, `genesets.gmt`, `run_summary.*`
- optional root `genesets.gmt` (combined)

Methods and equations: `docs/assays/drug_response/methods.tex`.

## Quickstart: generic long tables

```bash
geneset-extractors convert drug_response_screen \
  --response_tsv path/to/response.tsv \
  --drug_targets_tsv path/to/drug_targets.tsv \
  --out_dir results/drug_response \
  --organism human \
  --genome_build hg38

geneset-extractors validate results/drug_response
```

Minimal generic columns:

- response: `sample_id`, `drug_id`, `response`
- targets: `drug_id`, `gene_symbol` (optional `weight`, `source`)

Defaults:

- `response_metric=logfold_change`
- `response_direction=auto` (`logfold_change` -> `lower_is_more_sensitive`)
- `response_transform=robust_z_mad`
- `contrast_method=none` unless groups are available (then defaults to `group_vs_rest`)
- `scoring_model=target_weighted_sum`
- `ubiquity_penalty=fraction_active`
- `polypharm_downweight=true`
- `max_programs=50`
- `select=top_k`, `top_k=200`, `normalize=within_set_l1`
- strict GMT bounds: `gmt_min_genes=100`, `gmt_max_genes=500`

## Quickstart: PRISM convenience mode

```bash
geneset-extractors convert drug_response_screen \
  --prism_matrix_csv path/to/replicate_collapsed_logfold_change.csv \
  --prism_treatment_info_csv path/to/replicate_collapsed_treatment_info.csv \
  --prism_cell_line_info_csv path/to/cell_line_info.csv \
  --out_dir results/prism \
  --organism human \
  --genome_build hg38
```

PRISM defaults:

- sample_id from `row_name`
- drug_id from `broad_id` (fallback `column_name`)
- target parsing from treatment-info `target` column
- groups default to `primary_tissue` if `--groups_tsv` is not provided
- default contrast resolves to `group_vs_rest` when groups are available

## Contrast modes

- `none`: one program per sample (vulnerability profile)
- `group_mean`: one program per group (mean sensitivity)
- `group_vs_rest`: one signed program per group
- `case_control`: signed case-vs-control program(s)
  - requires `--case_control_tsv`
  - use `--case_control_within_group true` to compute per-group case-vs-control

When `gmt_split_signed` is not explicitly set, defaults to:

- `true` for `group_vs_rest` and `case_control`
- `false` otherwise

## Program count control

Screens can yield many programs. `--max_programs` defaults to 50.
If more are requested, converter truncates and warns:

- `group_vs_rest`/`case_control`/`group_mean`: rank by mean absolute drug contrast magnitude
- `none`: rank by number of non-missing sample-drug measurements

## Optional resources and overrides

Optional small resources:

- `--target_aliases_tsv` (or `--target_aliases_resource_id`)
- `--drug_blacklist_tsv` (or `--drug_blacklist_resource_id`)

Resource manager flags:

- `--resources_manifest`
- `--resources_dir`
- `--resource_policy {skip,fail}`

Bundled optional resource IDs:

- `drug_response_target_aliases_human`
- `drug_response_blacklist_human`

## Common warnings and fixes

- Low target coverage (`<50%` drugs mappable): provide better target annotations or alias mapping.
- Many invalid target tokens: inspect source target strings; provide alias map.
- Tiny groups in contrasts: increase group sample sizes or switch contrast mode.
- Program truncation by `--max_programs`: increase it for broader export.
- Empty/small GMT output: lower `--gmt_min_genes` or set `--emit_small_gene_sets true` for diagnostics.

## Output interpretation

- Positive signed sets (`__pos`) indicate relative sensitivity-associated target programs.
- Negative signed sets (`__neg`) indicate relative resistance-associated programs.
- For non-contrast/sample programs, scores are nonnegative vulnerability summaries by default.

## CLI to methods map

- `--response_metric`, `--response_direction`, `--response_transform` -> response orientation and normalization
- `--contrast_method`, `--case_control_within_group` -> program axis
- `--scoring_model` -> Model A / Model B
- `--ubiquity_penalty`, `--ubiquity_tau`, `--ubiquity_epsilon` -> ubiquity downweight
- `--polypharm_downweight`, `--polypharm_t0` -> polypharmacology downweight
- `--select`, `--top_k`, `--normalize`, `--emit_gmt`, GMT flags -> shared extraction/export
