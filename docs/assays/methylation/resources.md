# Methylation Resources

This page centralizes resource inputs for methylation converters.

## Required for probe-only CpG inputs

If `methylation_cpg_diff` input has only `probe_id` (no genomic coordinates),
provide one of:

1. `--probe_manifest_tsv <path>`
2. `--resources_dir` plus a resolvable default resource id based on
   `--array_type` and `--genome_build`

Expected resource filenames:

- `methylation_probe_manifest_450k_hg19.tsv.gz`
- `methylation_probe_manifest_epic_hg19.tsv.gz`
- `methylation_probe_manifest_450k_hg38.tsv.gz`
- `methylation_probe_manifest_epic_hg38.tsv.gz`

## Current manifest state

Bundled resource entries for methylation probe manifests are marked as
manual/bundle resources (no default URL in repo). This means `resources fetch`
will skip them unless URLs are supplied in an overlay manifest.

Canonical usage today is local bundle/direct mode:

1. Place manifest files under a local `--resources_dir`.
2. Set `OMICS2GENESET_RESOURCES_DIR` or pass `--resources_dir`.
3. Run converter without `--probe_manifest_tsv` and let auto-resolution select
   the build/array-specific id.

## Provenance checklist for shared bundles

When publishing or sharing local methylation manifest bundles, record:

- upstream source/package and exact version
- genome build
- transform steps used to create exported TSV
- license attribution used for redistributed derived files
