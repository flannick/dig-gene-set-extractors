# Proteomics PTM Reference Bundle

This guide describes the optional compact PTM bundle used by `ptm_site_diff` and `ptm_site_matrix`.

## Purpose

The v1 PTM bundle is intentionally lightweight. It does not try to merge raw quantitative phosphosite ratios across studies. It only provides:

1. site alias harmonization
2. a detection-frequency-based ubiquity prior

This is the safest public-bundle story for phosphoproteomics in the current repo.

## Bundle contents

Recommended human phosphoproteomics bundle files:

- `phosphosite_aliases_human_v1.tsv.gz`
- `phosphosite_ubiquity_human_v1.tsv.gz`
- `bundle_provenance.json`
- `local_resources_manifest.json`

The two runtime resource ids are:

- `phosphosite_aliases_human_v1`
- `phosphosite_ubiquity_human_v1`

Preset:

- `phosphoproteomics_default_optional_human`

## File schemas

### `phosphosite_aliases_human_v1.tsv.gz`

Required columns:

- `input_site_key`
- `canonical_site_key`
- `gene_id`
- `gene_symbol`
- `protein_accession`
- `residue`
- `position`
- `ptm_type`
- `source_dataset`

### `phosphosite_ubiquity_human_v1.tsv.gz`

Required columns:

- `canonical_site_key`
- `gene_id`
- `gene_symbol`
- `ptm_type`
- `n_samples_ref`
- `df_ref`
- `fraction_ref`
- `idf_ref`
- `n_datasets_ref`

## Build a local bundle

```bash
geneset-extractors workflows ptm_prepare_reference_bundle \
  --sources_tsv <sources.tsv> \
  --out_dir <ptm_bundle_dir> \
  --organism human \
  --ptm_type phospho \
  --bundle_id phosphoproteomics_human_v1
```

`<sources.tsv>` should contain at least:

- `path`
- `source_dataset`

Each source table should already be standardized enough that canonical site keys can be reconstructed from columns such as `site_id`, `protein_accession`, `residue`, and `position`.

## Runtime usage

### Direct bundle layout

If the bundle files are present directly in `--resources_dir`, the converter can resolve them without an overlay manifest:

```bash
geneset-extractors convert ptm_site_diff \
  --ptm_tsv <ptm_site_diff.tsv> \
  --out_dir <out_dir> \
  --organism human \
  --genome_build human \
  --resources_dir <ptm_bundle_dir>
```

The same bundle layout works for the matrix front-end:

```bash
geneset-extractors convert ptm_site_matrix \
  --ptm_matrix_tsv <ptm_matrix.tsv> \
  --sample_metadata_tsv <sample_metadata.tsv> \
  --out_dir <out_dir> \
  --organism human \
  --genome_build human \
  --resources_dir <ptm_bundle_dir>
```

### Local manifest overlay

The bundle workflow also writes `local_resources_manifest.json`. You can use it explicitly if you want the generated bundle to behave like a manifest overlay:

```bash
geneset-extractors convert ptm_site_diff \
  --ptm_tsv <ptm_site_diff.tsv> \
  --out_dir <out_dir> \
  --organism human \
  --genome_build human \
  --resources_dir <ptm_bundle_dir> \
  --resources_manifest <ptm_bundle_dir>/local_resources_manifest.json
```

## Missing-bundle behavior

Default behavior keeps the run easy:

- `--use_reference_bundle true`
- `--resource_policy skip`

If the alias or ubiquity file is missing, the converter:

- warns
- records the missing resource in metadata/run summaries
- falls back to direct canonicalization and `idf = 1`

For stricter runs:

```bash
--resource_policy fail
```

## Why the v1 bundle is compact

Public phosphoproteomics abundance tables are not safely comparable across studies as a shared atlas without additional within-study normalization and calibration. The v1 PTM bundle therefore avoids pooled raw abundance priors and keeps only:

- harmonization
- detection frequency / ubiquity weighting

That is intentional.
