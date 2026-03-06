# Gene Set Extractors Documentation Index

This is the repository-level documentation entrypoint.

Canonical Python package/module namespace is `geneset_extractors`.
CLI entrypoint aliases are provided for convenience.

## Namespaces and CLI names

- Python package namespace: `geneset_extractors`
- CLI aliases:
  - `geneset-extractors`
  - `geneset_extractors`

Both CLI names dispatch to the same implementation.

## Practical guides

- ATAC-seq: `docs/assays/atac/guide.md`
- RNA-seq: `docs/assays/rnaseq/guide.md`
- scRNA cNMF prep workflow: `docs/assays/rnaseq/scrna_cnmf_workflow.md`
- DNA methylation: `docs/assays/methylation/guide.md`
- DNA methylation resources: `docs/assays/methylation/resources.md`
- CNV segments: `docs/assays/cnv/guide.md`
- Drug response screens: `docs/assays/drug_response/guide.md`
- Drug response reference bundle: `docs/assays/drug_response/reference_bundle.md`
- Drug response compatibility guide: `docs/drug-response2geneset.md`

## Methods notes

- ATAC-seq methods: `docs/assays/atac/methods.tex`
- RNA-seq methods: `docs/assays/rnaseq/methods.tex`
- DNA methylation methods: `docs/assays/methylation/methods.tex`
- CNV methods: `docs/assays/cnv/methods.tex`
- Drug response methods: `docs/assays/drug_response/methods.tex`
- Methods index: `docs/methods.tex`

## Resource ID conventions and presets

Reference resources are managed via `src/geneset_extractors/resources/manifest.json` and optional overlay manifests.

- ID naming convention: `<assay>_<resource_kind>_<build>[_v<version>]`
  - Examples: `ccre_ubiquity_hg38`, `atac_reference_profiles_hg19`, `methylation_probe_manifest_450k_hg38`
- Presets group resource IDs by common workflows/builds (for example `atac_default_optional_hg19` and `atac_default_optional_hg38`).
- For resources without repository-hosted URLs, place files directly in `--resources_dir` using expected filenames, then run converters in direct/bundle mode.

Assay-specific resource instructions:

- ATAC resources and bundles: `docs/assays/atac/reference_bundle.md`
- Methylation probe manifests/resources: `docs/assays/methylation/resources.md`
- RNA and CNV are local-input-first (no required external bundle in current defaults).
- Drug response supports an optional local annotation/prior bundle for more connectable defaults.
- Drug response PRISM prep workflow: `geneset-extractors workflows prism_prepare` (see `docs/drug-response2geneset.md`)

Compatibility aliases remain at historical paths such as `docs/atac-seq2geneset.md` and `docs/rna-seq_methods.tex`.
