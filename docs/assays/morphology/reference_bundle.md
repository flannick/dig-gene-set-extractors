# Morphology Reference Bundles

Morphology extraction is most useful when the query is compared against a reusable perturbation reference library. This page describes the lightweight local bundle format used by `morphology_profile_query`.

## Bundle contents

A v1 bundle is a small directory containing text or parquet artifacts plus a manifest:

- `reference_profiles.tsv.gz` or `reference_profiles.parquet`
- `reference_metadata.tsv.gz`
- `compound_targets.tsv.gz`
- `feature_schema.tsv.gz`
- `feature_stats.tsv.gz`
- `<bundle_id>.bundle.json`

The bundle manifest maps logical roles to relative files:

- `reference_profiles`
- `reference_metadata`
- `compound_targets`
- `feature_schema`
- `feature_stats`

## How the converter uses a bundle

When you run:

```bash
geneset-extractors convert morphology_profile_query \
  --query_profiles_tsv path/to/query_profiles.tsv \
  --reference_bundle_id morphology_jump_target_pilot_u2os_48h_v1 \
  --resources_dir /path/to/bundle_root \
  --out_dir results/morph_query \
  --organism human \
  --genome_build hg38
```

the converter resolves `morphology_jump_target_pilot_u2os_48h_v1.bundle.json` from `--resources_dir`, then loads the referenced files relative to that manifest.

No separate resource-manifest overlay is required when the bundle manifest is present directly in `--resources_dir`.
Resolution order is:

1. local `--resources_dir/<bundle_id>.bundle.json`
2. resource-manager manifest lookup

## Resource ids

Current built-in manual/local bundle resource ids:

- `morphology_jump_target_pilot_u2os_48h_v1`
- `morphology_jump_target_pilot_a549_48h_v1`

These entries are defined in `src/geneset_extractors/resources/manifest.json` as local/manual resources. If you do not have a hosted URL, place the bundle files directly in `--resources_dir` with the expected manifest filename.

## Recommended workflow

1. Build or obtain a bundle directory.
2. Point `--resources_dir` at that directory.
3. Use `--reference_bundle_id` in the converter.
4. Keep bundle cell line and timepoint matched to the query experiment when possible.

## Building a bundle locally

Use the helper workflow:

```bash
geneset-extractors workflows jump_prepare_reference_bundle \
  --profile_paths path/to/profiles.tsv \
  --experimental_metadata_tsv path/to/experimental_metadata.tsv \
  --compound_targets_tsv path/to/jump_target_compound_metadata.tsv \
  --cell_type_filter U2OS \
  --timepoint_filter 48 \
  --out_dir results/jump_u2os_48h_bundle
```

This creates a local bundle manifest plus the referenced files.

The workflow accepts either:

- an already standardized `compound_targets.tsv`, or
- common public JUMP target metadata TSVs, which it normalizes internally to `compound_id,gene_symbol,weight,source`.

Bundle outputs now also include:

- `bundle_summary.json`
- `bundle_summary.txt`

These summarize which modalities were included, which were omitted, and why.

## Missing-bundle behavior

- `--resource_policy skip`:
  - warn and continue only if explicit reference files were supplied as fallback.
- `--resource_policy fail`:
  - missing bundle is a hard error.

Bundle usage and missing-resource information are recorded in `geneset.meta.json` and run summaries.

## Context coherence defaults

The bundle workflow defaults to same-timepoint coherence:

- `--require_same_timepoint_across_modalities true`
- `--allow_mixed_timepoints false`
- `--allow_missing_modalities true`

This means the builder will keep a coherent cell-type/timepoint context and warn if some modalities are unavailable, rather than silently mixing 48h and 96h references.

## Hubness penalty support

The bundle builder computes a `hub_score` for each reference perturbation from within-bundle reference-reference similarities. The extractor can use this score to downweight generic hub-like references:

- `--hubness_penalty inverse_rank` (default)
- `--hubness_penalty inverse_linear`
- `--hubness_penalty none`

The extractor also keeps only a narrow local neighborhood by default:

- `--max_reference_neighbors 20`
- `--min_specificity_confidence_to_emit_opposite medium`

This is deliberate. Morphology retrieval becomes much less specific when weak long-tail neighbors are allowed to contribute.
