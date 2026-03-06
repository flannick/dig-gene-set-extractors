# Drug Response Reference Bundle

`drug_response_screen` can run without any external bundle if you already have:

- a response table or PRISM inputs, and
- a usable drug-to-target mapping.

A local drug-response reference bundle improves defaults for PRISM/GDSC/CTRP-like screens by
adding stable target annotations, drug ID normalization, target ubiquity priors, and compound QC.

## Bundle contents

Recommended bundle files for the current `connectable` preset:

- `drug_target_edges_human_v1.tsv.gz`
  - columns: `drug_id`, `gene_symbol`, `weight`, `confidence`, `source`
- `drug_alias_map_human_v1.tsv.gz`
  - columns: `input_drug_id`, `canonical_drug_id`, `namespace`, `source`
- `target_ubiquity_human_v1.tsv.gz`
  - columns: `gene_symbol`, `n_drugs`, `idf`, optional `family`
- `compound_qc_human_v1.tsv.gz`
  - columns: `drug_id`, `n_targets`, `polypharm_flag`, `pan_toxic_flag`, `recommended_use`, `blacklist_reason`

These resource ids are registered in `src/geneset_extractors/resources/manifest.json`:

- `drug_target_edges_human_v1`
- `drug_alias_map_human_v1`
- `target_ubiquity_human_v1`
- `compound_qc_human_v1`

Current status:

- the repo supports local bundle resolution from `--resources_dir`
- the bundled manifest contains placeholder metadata for these files
- if you have a local bundle, place the files at the bundle root using the exact filenames above

## How to use a local bundle

Point `--resources_dir` at the directory containing the bundle files:

```bash
geneset-extractors convert drug_response_screen \
  --response_tsv path/to/response.tsv \
  --out_dir results/drug_response \
  --organism human \
  --genome_build hg38 \
  --resources_dir /path/to/drug_response_bundle
```

Or export the shared environment variable:

```bash
export GENESET_EXTRACTORS_RESOURCES_DIR=/path/to/drug_response_bundle
```

Then run the converter without repeating `--resources_dir`.

## What changes when the bundle is present

Default preset is `--program_preset connectable`.
When bundle resources are present, the converter will try to:

1. normalize drug IDs with `drug_alias_map_human_v1`
2. use `drug_target_edges_human_v1` if `--drug_targets_tsv` is not provided
3. use `target_ubiquity_human_v1` for a stable target-side IDF penalty
4. use `compound_qc_human_v1` to:
   - drop compounds with `recommended_use=drop` or `pan_toxic_flag=true`
   - warn on compounds with `recommended_use=warn`
   - downweight compounds flagged `polypharm_flag=true`

Without the bundle, the extractor still runs, but defaults are broader:

- target ubiquity is derived from the current retained compound subset
- PRISM mode falls back to treatment-info target text when no target-edge table is available
- nuisance compounds are not bundle-filtered

## Fallback behavior and warnings

Under `--resource_policy skip` (default), missing bundle files do not stop the run.
Instead the converter records them in metadata and warns with messages like:

- `warning: bundle_resource_missing resource=target_ubiquity_human_v1 action=compute_subset_derived_target_ubiquity`
- `warning: bundle_resource_missing resource=drug_target_edges_human_v1 action=fall_back_to_prism_treatment_info_or_error`

Under `--resource_policy fail`, missing requested bundle files raise an error.

## Metadata and reproducibility

When bundle files are used, they are recorded in:

- `program=*/geneset.meta.json`
- root `run_summary.json`

Look under:

- `summary.target_summary.resources.used`
- `summary.target_summary.resources.missing`

This lets you distinguish:

- explicit user-provided target annotations
- bundle-derived target edges / ubiquity priors / compound QC
- fallback to PRISM treatment-info targets

## Recommended defaults vs exploratory mode

- Recommended: `--program_preset connectable`
  - uses bundle priors when available
  - aims for more connectable, less library-driven outputs
- Exploratory: `--program_preset broad_pharmacology`
  - lighter-touch filtering and warnings
  - useful when broad pharmacology itself is the object of interest
