# Calorimetry Reference Bundles

`calr_profile_query` is most useful when a query can be compared against a reusable calorimetry reference library. `workflows calr_prepare_reference_bundle` builds that local bundle in the same pattern used by the morphology bundle workflow: one unpacked bundle directory plus an optional distribution tarball.

## Bundle contents

A calorimetry bundle contains:

- `reference_profiles.tsv.gz`
- `reference_metadata.tsv.gz`
- `feature_schema.tsv.gz`
- `feature_stats.tsv.gz`
- `term_templates.tsv.gz`
- `phenotype_gene_edges.tsv.gz`
- optional `term_hierarchy.tsv.gz`
- `<bundle_id>.bundle.json`
- `bundle_summary.json`
- `bundle_summary.txt`

The bundle manifest maps logical roles to these relative files.

Profile-query uses:

- `reference_profiles`
- `reference_metadata`
- optional `feature_schema`
- optional `feature_stats`

Ontology-mapper can also use the same bundle for:

- `term_templates`
- `phenotype_gene_edges`
- optional `term_hierarchy`

## Packaged defaults vs bundle resources

`calr_ontology_mapper` can run without a bundle because mouse packaged defaults are included in the package:

- `calr_term_templates_mouse_v1.tsv`
- `calr_phenotype_gene_edges_mouse_v1.tsv`
- `calr_term_hierarchy_mouse_v1.tsv`

`calr_profile_query` cannot run without either:

- explicit reference tables, or
- a bundle containing the reference profiles and metadata

For human runs, do not rely on the packaged ontology defaults. They are mouse-first resources. Provide explicit human resources or a human bundle.

## Build a bundle locally

```bash
geneset-extractors workflows calr_prepare_reference_bundle \
  --reference_profiles_tsv <reference_profiles.tsv> \
  --reference_metadata_tsv <reference_metadata.tsv> \
  --out_dir <bundle_dir> \
  --organism mouse \
  --bundle_id calorimetry_mouse_v1
```

Optional explicit inputs:

- `--feature_schema_tsv`
- `--feature_stats_tsv`
- `--term_templates_tsv`
- `--phenotype_gene_edges_tsv`
- `--term_hierarchy_tsv`

If the schema/stats are not provided, the workflow derives them from `reference_profiles_tsv`.
If the ontology files are not provided, the workflow falls back to the packaged mouse defaults.

## Distribution artifact

By default the workflow also writes:

- `dist/dig-gene-set-extractors-<bundle_id>.tar.gz`
- `dist/SHA256SUMS.txt`
- `dist/distribution_manifest.json`

This supports the same user story as the morphology bundle flow:

1. build or download one tarball
2. extract it
3. point `--resources_dir` at the extracted bundle directory

## Using a bundle

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

Bundle resolution follows the local-manifest pattern used elsewhere in the repo:

1. `<resources_dir>/<bundle_id>.bundle.json`
2. resource-manifest lookup if a manifest overlay is provided

## Metadata expected in the reference library

The bundle is more useful when `reference_metadata.tsv` records at least:

- primary perturbed gene or reference id
- perturbation type
- selected mass covariate
- acclimation state
- ambient temperature
- sex
- site or run identifier
- QC or hubness-related fields when available

These metadata fields are used for:

- provenance mismatch penalties
- confidence summaries
- bundle interpretation

## Missing-bundle behavior

- `--resource_policy skip`: warn and continue only if enough explicit fallback files are supplied
- `--resource_policy fail`: missing required bundle files are a hard error
