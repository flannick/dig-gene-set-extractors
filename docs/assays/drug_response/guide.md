# Drug Response Screens to Gene Sets (`drug_response_screen`)

This guide covers conversion of drug response screens (PRISM, GDSC/CTRP-like long tables)
into gene-weighted programs and GMT gene sets.

`drug_response_screen` builds gene programs from:

- response values per sample and drug
- drug-target annotations

and emits grouped outputs:

- `manifest.tsv`
- `group_qc.tsv` (all groups including skipped with reasons)
- root `run_summary.json` / `run_summary.txt`
- `program=<PROGRAM>/geneset.tsv`
- `program=<PROGRAM>/geneset.meta.json`
- optional per-program `geneset.full.tsv`, `genesets.gmt`, `run_summary.*`
- optional root `genesets.gmt` (combined)

Methods and equations: `docs/assays/drug_response/methods.tex`.
Optional local bundle: `docs/assays/drug_response/reference_bundle.md`.

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

- `program_preset=connectable`
- `response_metric=logfold_change`
- `response_direction=auto` (`logfold_change` -> `lower_is_more_sensitive`)
- `response_transform=robust_z_mad`
- `contrast_method=none` unless groups are available (then defaults to `group_vs_rest`)
- `scoring_model=target_weighted_sum`
- `response_ubiquity_penalty=fraction_active` (legacy alias: `ubiquity_penalty`)
- `target_ubiquity_penalty=idf`
- `polypharm_downweight=true`
- `min_group_size=5`
- `max_targets_per_drug=50`, `target_promiscuity_policy=warn`
- `max_programs=50`
- `select=top_k`, `top_k=200`, `normalize=within_set_l1`
- strict GMT bounds: `gmt_min_genes=100`, `gmt_max_genes=500`

If a local bundle is available via `--resources_dir` (or `GENESET_EXTRACTORS_RESOURCES_DIR`),
`connectable` will also try to:

- normalize drug IDs through `drug_alias_map_human_v1`
- use `drug_target_edges_human_v1` when `--drug_targets_tsv` is omitted
- use `target_ubiquity_human_v1` instead of deriving target IDF from the current subset
- apply `compound_qc_human_v1` nuisance-compound filtering/warnings

## Quickstart: PRISM convenience mode

Recommended first step (download + standardize into long tables):

```bash
geneset-extractors workflows prism_prepare \
  --out_dir results/prism_prepare
```

Default file-id template is Figshare (`https://ndownloader.figshare.com/files/{file_id}`).
`prism_prepare` performs content sniffing on downloads and warns/fails if content
looks like portal metadata JSON/HTML instead of tabular PRISM data.

Quick deterministic subset for fast iteration:

```bash
geneset-extractors workflows prism_prepare \
  --out_dir results/prism_prepare_quick \
  --subset_seed 7 \
  --balance_by primary_tissue \
  --min_per_balance_bin 5 \
  --max_cell_lines_per_group 50 \
  --max_cell_lines_total 500 \
  --max_compounds_total 300
```

Then run the converter on standardized tables:

```bash
geneset-extractors convert drug_response_screen \
  --response_tsv results/prism_prepare/response_long.tsv \
  --drug_targets_tsv results/prism_prepare/drug_targets.tsv \
  --groups_tsv results/prism_prepare/groups.tsv \
  --out_dir results/prism \
  --organism human \
  --genome_build hg38
```

With a local bundle for more robust defaults:

```bash
geneset-extractors convert drug_response_screen \
  --response_tsv results/prism_prepare/response_long.tsv \
  --groups_tsv results/prism_prepare/groups.tsv \
  --out_dir results/prism_connectable \
  --organism human \
  --genome_build hg38 \
  --resources_dir /path/to/drug_response_bundle
```

Direct PRISM input mode is also supported:

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

## Group and target QC defaults

- `--min_group_size 5`: skip tiny groups by default with explicit warnings and metadata.
- `--max_targets_per_drug 50` and `--target_promiscuity_policy {warn,drop,cap}`:
  control very large target lists.
- `--target_ubiquity_penalty idf`:
  downweights genes targeted by many retained drugs.

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

Drug-response reference bundle IDs:

- `drug_target_edges_human_v1`
- `drug_alias_map_human_v1`
- `target_ubiquity_human_v1`
- `compound_qc_human_v1`

See `docs/assays/drug_response/reference_bundle.md` for the expected bundle layout and fallback behavior.

## Common warnings and fixes

- Low target coverage (`<50%` drugs mappable): provide better target annotations or alias mapping.
- Many invalid target tokens: inspect source target strings; provide alias map.
- Tiny groups in contrasts: increase group sample sizes or switch contrast mode.
- Tiny groups skipped by default: lower `--min_group_size` only if biologically justified.
- Program truncation by `--max_programs`: increase it for broader export.
- Empty/small GMT output: lower `--gmt_min_genes` or set `--emit_small_gene_sets true` for diagnostics.
- Promiscuous targets: adjust `--max_targets_per_drug` and `--target_promiscuity_policy`.
- Broad GPCR/neuroactive signal: keep `--target_ubiquity_penalty idf` and inspect run summaries.
- Without the bundle, outputs can be broader and more library-driven because target ubiquity and nuisance compound handling are subset-derived or absent.
- Downloaded JSON/HTML instead of PRISM tables: inspect `prepare_summary.json` `fetch.*.sniff`;
  use Figshare template `https://ndownloader.figshare.com/files/{file_id}`.

## GMT format

- Default: `--gmt_format dig2col` (`set_id<TAB>gene1 gene2 ...`)
- Classic parser compatibility: `--gmt_format classic` (`set_id<TAB>na<TAB>gene1<TAB>gene2...`)

## Output interpretation

- Positive signed sets (`__pos`) indicate relative sensitivity-associated target programs.
- Negative signed sets (`__neg`) indicate relative resistance-associated programs.
- For non-contrast/sample programs, scores are nonnegative vulnerability summaries by default.

## CLI to methods map

- `--response_metric`, `--response_direction`, `--response_transform` -> response orientation and normalization
- `--contrast_method`, `--case_control_within_group` -> program axis
- `--scoring_model` -> Model A / Model B
- `--response_ubiquity_penalty` / `--ubiquity_penalty`, `--ubiquity_tau`, `--ubiquity_epsilon` -> response-side ubiquity downweight
- `--target_ubiquity_penalty` -> target-side ubiquity downweight
- `--polypharm_downweight`, `--polypharm_t0` -> polypharmacology downweight
- `--max_targets_per_drug`, `--target_promiscuity_policy` -> promiscuity handling
- `--select`, `--top_k`, `--normalize`, `--emit_gmt`, GMT flags -> shared extraction/export
