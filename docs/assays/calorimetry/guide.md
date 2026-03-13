# Calorimetry Guide

This assay family is for CalR-compatible indirect calorimetry datasets. It is a phenotype-to-gene extractor family, not a direct molecular assay. The extractor first summarizes physiology in a CalR-aware way, then routes that physiology to genes through either phenotype terms or a gene-labeled reference library.

The calorimetry physiology and packaged ontology resources are mouse-first, but the default emitted gene lists are now humanized. By default:

- source mouse genes are mapped to human orthologs
- only unique mouse-to-human mappings are kept
- source mouse gene provenance remains in metadata and public-study reference tables

Current public entrypoints:

- `calr_ontology_mapper`
- `calr_profile_query`
- `workflows calr_prepare_public`
- `workflows calr_prepare_reference_bundle`

The implementation is mouse-first. The packaged ontology defaults are mouse phenotype resources. For human runs, provide explicit ontology resources or a human-calibrated reference bundle instead of relying on packaged defaults.

## Quickstart

### Ontology mapper with packaged mouse defaults

```bash
geneset-extractors convert calr_ontology_mapper \
  --calr_data_csv <calr_data.csv> \
  --session_csv <session.csv> \
  --out_dir <out_dir> \
  --organism mouse \
  --genome_build mm39
```

This uses packaged default files for:

- term templates
- phenotype-gene edges
- term hierarchy
- mouse-to-human ortholog table used to humanize emitted genes

### Profile query with explicit reference tables

```bash
geneset-extractors convert calr_profile_query \
  --calr_data_csv <calr_data.csv> \
  --session_csv <session.csv> \
  --reference_profiles_tsv <reference_profiles.tsv> \
  --reference_metadata_tsv <reference_metadata.tsv> \
  --feature_schema_tsv <feature_schema.tsv> \
  --feature_stats_tsv <feature_stats.tsv> \
  --out_dir <out_dir> \
  --organism mouse \
  --genome_build mm39
```

By default this emits human gene programs even when the reference metadata is still mouse-first. Override with:

- `--output_gene_species source`
- `--ortholog_policy expand_all`
- `--mouse_human_orthologs_tsv <custom.tsv.gz>`

### Build and use a local calorimetry reference bundle

Build:

```bash
geneset-extractors workflows calr_prepare_reference_bundle \
  --reference_profiles_tsv <reference_profiles.tsv> \
  --reference_metadata_tsv <reference_metadata.tsv> \
  --out_dir <bundle_dir> \
  --organism mouse \
  --bundle_id calorimetry_mouse_v1
```

Query:

```bash
geneset-extractors convert calr_profile_query \
  --calr_data_csv <calr_data.csv> \
  --session_csv <session.csv> \
  --reference_bundle_id calorimetry_mouse_v1 \
  --resources_dir <bundle_dir> \
  --out_dir <out_dir> \
  --organism mouse \
  --genome_build mm39
```

`calr_ontology_mapper` can also use a bundle when the bundle includes:

- `term_templates`
- `phenotype_gene_edges`
- optional `term_hierarchy`

### Build a profile-query bundle from raw public Cal-Repository studies

First create a small `studies.tsv` manifest with one row per reference study/group. Required columns are:

- `study_id`
- `calr_data_csv`
- `session_csv`
- `reference_group`
- `gene_symbol`

Minimal example:

```tsv
study_id\tcalr_data_csv\tsession_csv\treference_group\tgene_symbol
FIG2_HFD\t/path/Fig2_CalR_data.csv\t/path/Fig2_CalR_Davis_Session.csv\tHFD\tLepr
```

Then run:

```bash
geneset-extractors workflows calr_prepare_public \
  --studies_tsv <studies.tsv> \
  --out_dir <public_bundle_out> \
  --organism mouse \
  --bundle_id calorimetry_public_mouse_v1
```

This writes:

- `reference_profiles.tsv`
- `reference_metadata.tsv`
- `feature_schema.tsv`
- `feature_stats.tsv`
- `resolved_studies.tsv`
- `prepare_summary.json`
- `bundle/<bundle_id>.bundle.json`

and, by default for mouse studies, humanized reference metadata fields:

- `output_gene_id`
- `output_gene_symbol`
- `output_gene_species`
- `source_gene_id`
- `source_gene_symbol`
- `source_gene_species`
- `gene_mapping_status`

You can then query that bundle directly with `calr_profile_query`.

## Input expectations

The public converters treat a standardized CalR data CSV as the native input.

Preferred inputs:

- `--calr_data_csv`
- `--session_csv`

Optional:

- `--exclusions_tsv`

Supported session layouts now include both:

- direct row-per-subject session tables with explicit group/window columns
- authentic wide CalR session layouts where `group1`, `group2`, ... contain subject membership, `group_names` provides the display labels, and `xrange` / `light` / `exc` define the analysis window and photoperiod

For wide CalR layouts the extractor derives and records:

- `session_group_layout_mode = wide_membership`
- `session_window_layout_mode = calr_matrix`

For mouse-first ontology or reference bundles, the default output routing also records:

- `output_gene_species = human`
- `orthology_summary`

The session file is strongly preferred. If it is absent, the extractor can still run only in exploratory mode. That mode emits a hard warning and records:

- `session_mode = exploratory`

If a session file is present, the output records:

- `session_mode = explicit`
- `session_group_layout_mode`
- `session_window_layout_mode`

## CalR processing rules preserved by the extractor

These are the important semantics that the implementation keeps from CalR:

- selected analysis window is not silently the full trace
- session-provided window and photoperiod information take priority over fallback heuristics
- no naive body-weight normalization as the primary analysis
- mass-dependent variables use explicit mass/body-composition covariate logic
- mass-independent variables do not use that covariate by default
- `ee` and `rer` are used as reported when present; they are not silently recomputed
- interaction with the mass covariate is tested and retained only when warranted
- exclusions are explicit and auditable
- authentic CalR session matrices are not treated as subject-vs-rest group labels

Current covariate preference order:

1. explicit `--mass_covariate`
2. `Lean.Mass`
3. `Total.Mass`
4. `subject.mass`
5. no covariate only as an exploratory fallback with warning

## Program outputs

Both calorimetry converters emit grouped outputs across contrasts and physiology programs.

Current program decomposition:

- `global`
- `thermogenesis`
- `substrate_use`
- `intake_balance`
- `activity_circadian`

`calr_ontology_mapper` emits:

- `mode=core`
- `mode=expanded`

`calr_profile_query` currently emits:

- `mode=core`

Each emitted program directory contains:

- `geneset.tsv`
- `geneset.meta.json`
- optional `geneset.full.tsv`
- optional `genesets.gmt`
- `run_summary.json`
- `run_summary.txt`

At the dataset root:

- `manifest.tsv`
- optional combined `genesets.gmt`
- root `run_summary.json`
- root `run_summary.txt`

## Reference retrieval behavior

`calr_profile_query` aligns query features against a reference schema and standardizes with bundle or explicit feature stats. It then scores references with:

- `--similarity_metric cosine|pearson`
- `--similarity_floor`
- `--similarity_power`
- `--hubness_penalty none|inverse_linear|inverse_rank`
- `--provenance_mismatch_penalty`

Low-confidence retrieval is warned explicitly. Provenance mismatches such as different mass-covariate or acclimation metadata are penalized, not silently ignored.

When a reference bundle already carries humanized `output_gene_*` fields, those are used directly. Otherwise the profile-query output is humanized after retrieval using the bundled or packaged mouse-human ortholog table.

## Selection and GMT defaults

Calorimetry follows the repo-wide conservative GMT defaults:

- `--select top_k`
- `--top_k 200`
- `--normalize within_set_l1`
- `--emit_gmt true`
- `--gmt_topk_list 200`
- `--gmt_min_genes 100`
- `--gmt_max_genes 500`
- `--emit_small_gene_sets false`

## Common pitfalls

- No session file: the run is exploratory and may not match a CalR analysis window.
- Human runs with packaged defaults: packaged ontology resources are mouse-only; use explicit human resources or a human bundle.
- Ambiguous orthologs: default `--ortholog_policy unique_only` drops ambiguous mouse-to-human mappings rather than guessing.
- Full-trace analysis: if acclimation is included, the output may not match a typical post-acclimation calorimetry analysis.
- Missing body-composition covariate: mass-dependent variables will fall back and warn.
- Combined or crossover designs: v1 is intended for straightforward baseline or chronic grouped designs, not silent crossover handling.

## When to use each converter

Use `calr_ontology_mapper` when:

- you want the lowest-friction phenotype-to-gene path
- you do not yet have a curated reference library
- you want core and expanded phenotype-routed outputs

Use `calr_profile_query` when:

- you have a local calorimetry reference library
- you want nearest-profile evidence with provenance-aware penalties
- you want to reuse the same reference bundle across many studies

Use `workflows calr_prepare_public` when:

- your starting point is raw local Cal-Repository-style `*_data.csv` plus `*_Session.csv` files
- you want a reproducible bridge from public studies to a gene-labeled `calr_profile_query` bundle
- you want the bundle to preserve CalR window/group semantics rather than hand-building reference profiles
