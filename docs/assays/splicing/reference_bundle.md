# Alternative Splicing Reference Bundle

This guide describes the optional compact splicing bundle used by `splice_event_diff` and `splice_event_matrix`.

## Purpose

The v1 splicing bundle is intentionally compact. It does not try to pool raw PSI distributions across unrelated studies for direct biological scoring. It only provides:

1. event alias harmonization
2. a detection-frequency ubiquity prior
3. a conservative impact prior that shrinks toward neutral
4. an optional per-gene event-burden reference count

## Bundle contents

Recommended human splicing bundle files:

- `splice_event_aliases_human_v1.tsv.gz`
- `splice_event_ubiquity_human_v1.tsv.gz`
- `splice_event_impact_human_v1.tsv.gz`
- `splice_gene_event_burden_human_v1.tsv.gz`
- `bundle_provenance.json`
- `local_resources_manifest.json`

Runtime resource ids:

- `splice_event_aliases_human_v1`
- `splice_event_ubiquity_human_v1`
- `splice_event_impact_human_v1`
- `splice_gene_event_burden_human_v1`

Preset:

- `splicing_default_optional_human`

## File schemas

### `splice_event_aliases_human_v1.tsv.gz`

Required columns:

- `input_event_key`
- `canonical_event_key`
- `canonicalization_status`
- `canonicalization_confidence`
- `gene_id`
- `gene_symbol`
- `event_type`
- `chrom`
- `start`
- `end`
- `strand`
- `source_dataset`

### `splice_event_ubiquity_human_v1.tsv.gz`

Required columns:

- `canonical_event_key`
- `gene_id`
- `gene_symbol`
- `event_type`
- `canonicalization_status`
- `canonicalization_confidence`
- `n_samples_ref`
- `df_ref`
- `fraction_ref`
- `idf_ref`
- `n_datasets_ref`

### `splice_event_impact_human_v1.tsv.gz`

Required columns:

- `canonical_event_key`
- `gene_id`
- `gene_symbol`
- `event_type`
- `canonicalization_status`
- `canonicalization_confidence`
- `impact_weight_raw`
- `impact_evidence`
- `annotation_status`

### `splice_gene_event_burden_human_v1.tsv.gz`

Required columns:

- `gene_symbol`
- `n_canonical_events_ref`
- `n_high_confidence_events_ref`
- `n_low_confidence_events_ref`

## Build a local bundle

```bash
geneset-extractors workflows splice_prepare_reference_bundle \
  --sources_tsv <sources.tsv> \
  --out_dir <splice_bundle_dir> \
  --organism human \
  --bundle_id splice_human_v1
```

`<sources.tsv>` must contain at least:

- `path`
- `source_dataset`

Each `path` should point to a staged `bundle_source_row.tsv`, usually produced by:

```bash
geneset-extractors workflows splice_prepare_public \
  --input_mode tcga_spliceseq \
  --psi_tsv <tcga_spliceseq.tsv> \
  --sample_annotations_tsv <sample_annotations.tsv> \
  --out_dir <prepared_dir> \
  --organism human \
  --genome_build hg38 \
  --study_id <study_id>
```

## Runtime usage

If the bundle files are present directly in `--resources_dir`, the converters can resolve them without an overlay manifest:

```bash
geneset-extractors convert splice_event_diff \
  --splice_tsv <splice_event_diff.tsv> \
  --out_dir <out_dir> \
  --organism human \
  --genome_build hg38 \
  --resources_dir <splice_bundle_dir>
```

The same bundle layout works for the matrix front-end:

```bash
geneset-extractors convert splice_event_matrix \
  --psi_matrix_tsv <psi_matrix.tsv> \
  --sample_metadata_tsv <sample_metadata.tsv> \
  --out_dir <out_dir> \
  --organism human \
  --genome_build hg38 \
  --resources_dir <splice_bundle_dir>
```

If `--resources_dir` contains `local_resources_manifest.json`, the splicing converters auto-discover and merge that local manifest even when `--resources_manifest` is not passed explicitly.

## Canonicalization confidence and low-confidence keys

The bundle now carries explicit canonicalization confidence:

- `coordinate_canonical` / `high`
- `raw_id_fallback` / `low`

Interpretation rule:

- high-confidence coordinate keys are the closest thing this assay family has to cross-study harmonized event identities
- low-confidence raw-ID keys are cohort-local and should not be interpreted as globally harmonized events

Runtime implications:

- alias resolution still allows exact-string matches on low-confidence raw ids
- broader cross-study collapsing is intentionally avoided for those keys
- ubiquity and impact priors are shrunk strongly toward neutral when only low-confidence matching is available

## Missing-bundle behavior

Default runtime behavior keeps the extractor easy to use:

- `--use_reference_bundle true`
- `--resource_policy skip`

If alias, ubiquity, or impact resources are missing, the converter:

- warns
- records the missing resource in metadata and run summaries
- falls back to direct canonicalization and neutral priors

For stricter runs:

```bash
--resource_policy fail
```
