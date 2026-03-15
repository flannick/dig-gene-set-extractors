# Alternative Splicing Reference Bundle

This guide describes the optional compact splicing bundle used by `splice_event_diff` and `splice_event_matrix`.

## Purpose

The v1 splicing bundle is intentionally compact. It does not try to pool raw PSI distributions across unrelated studies for direct biological scoring. It only provides:

1. event alias harmonization
2. a detection-frequency ubiquity prior
3. a conservative impact prior that shrinks toward neutral
4. an optional per-gene event-burden reference count
5. dataset-stratified nuisance summaries used to exclude same-cohort priors at runtime
6. an optional compact gene-locus table for chromosome/arm-aware locus suppression

## Bundle contents

Recommended human splicing bundle files:

- `splice_event_aliases_human_v1.tsv.gz`
- `splice_event_ubiquity_human_v1.tsv.gz`
- `splice_event_ubiquity_by_dataset_human_v1.tsv.gz`
- `splice_event_impact_human_v1.tsv.gz`
- `splice_gene_event_burden_human_v1.tsv.gz`
- `splice_gene_event_burden_by_dataset_human_v1.tsv.gz`
- optional `splice_gene_locus_human_v1.tsv.gz`
- `bundle_provenance.json`
- `local_resources_manifest.json`

Runtime resource ids:

- `splice_event_aliases_human_v1`
- `splice_event_ubiquity_human_v1`
- `splice_event_ubiquity_by_dataset_human_v1`
- `splice_event_impact_human_v1`
- `splice_gene_event_burden_human_v1`
- `splice_gene_event_burden_by_dataset_human_v1`
- `splice_gene_locus_human_v1`

Preset:

- `splicing_default_optional_human`

## File schemas

### `splice_event_aliases_human_v1.tsv.gz`

Required columns:

- `input_event_key`
- `canonical_event_key`
- `canonicalization_status`
- `canonicalization_confidence`
- `event_key_namespace`
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
- `event_key_namespace`
- `n_samples_ref`
- `df_ref`
- `fraction_ref`
- `idf_ref`
- `n_datasets_ref`
- `fraction_datasets_ref`

### `splice_event_ubiquity_by_dataset_human_v1.tsv.gz`

Required columns:

- `canonical_event_key`
- `source_dataset`
- `gene_id`
- `gene_symbol`
- `event_type`
- `canonicalization_status`
- `canonicalization_confidence`
- `event_key_namespace`
- `n_samples_ref`
- `df_ref`

### `splice_event_impact_human_v1.tsv.gz`

Required columns:

- `canonical_event_key`
- `gene_id`
- `gene_symbol`
- `event_type`
- `canonicalization_status`
- `canonicalization_confidence`
- `event_key_namespace`
- `impact_weight_raw`
- `impact_evidence`
- `annotation_status`
- `n_datasets_ref`

### `splice_gene_event_burden_human_v1.tsv.gz`

Required columns:

- `gene_symbol`
- `n_canonical_events_ref`
- `n_high_confidence_events_ref`
- `n_medium_confidence_events_ref`
- `n_low_confidence_events_ref`
- `n_unique_event_groups_ref`
- `n_studies_ref`
- `n_studies_high_confidence_ref`
- `fraction_low_confidence_events_ref`
- `median_unique_groups_per_study`

### `splice_gene_event_burden_by_dataset_human_v1.tsv.gz`

Required columns:

- `gene_symbol`
- `source_dataset`
- `n_canonical_events_ref`
- `n_high_confidence_events_ref`
- `n_medium_confidence_events_ref`
- `n_low_confidence_events_ref`
- `n_unique_event_groups_ref`

### Optional `splice_gene_locus_human_v1.tsv.gz`

Recommended columns:

- `gene_symbol`
- `gene_id`
- `chrom`
- `start`
- `end`
- `chromosome_arm`

## Build a local bundle

```bash
geneset-extractors workflows splice_prepare_reference_bundle \
  --sources_tsv <sources.tsv> \
  --out_dir <splice_bundle_dir> \
  --organism human \
  --bundle_id splice_human_v1
```

For leave-one-cohort-out validation, exclude the target cohort when building the prior bundle:

```bash
geneset-extractors workflows splice_prepare_reference_bundle \
  --sources_tsv <sources.tsv> \
  --out_dir <splice_bundle_dir> \
  --organism human \
  --bundle_id splice_human_v1 \
  --exclude_source_datasets TCGA_BRCA
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

## Confidence tiers and namespaces

The bundle now carries explicit canonicalization confidence:

- `coordinate_canonical` / `high`
- `source_family_stable_id` / `medium`
- `raw_id_fallback` / `low`

And an explicit event namespace:

- `global_coordinate` for coordinate-canonical keys
- `tcga_spliceseq_asid` for TCGA SpliceSeq source-family-stable AS identifiers
- `raw_id_fallback` for generic cohort-local fallbacks

Interpretation rule:

- high-confidence coordinate keys are the closest thing this assay family has to cross-study harmonized event identities
- medium-confidence `tcga_spliceseq_asid` keys are reusable across TCGA SpliceSeq-like datasets only
- low-confidence raw-ID keys are cohort-local and should not be interpreted as globally harmonized events

Runtime implications:

- alias resolution still allows exact-string matches on low-confidence raw ids
- medium-confidence priors are only used within the matching namespace
- broader cross-study collapsing is intentionally avoided for low-confidence keys
- low-confidence ubiquity priors are neutralized to exactly `1.0`
- low-confidence impact priors are only allowed to modulate scores when the bundle reports recurrence across at least two source datasets

The same raw event ID in a generic diff table is therefore not promoted to `medium`; the source-family-stable tier is only for TCGA SpliceSeq-like inputs.

## Study-level nuisance summaries

The bundle now treats study identity as part of the nuisance model rather than only pooled sample opportunity.

Current v1 summaries include:

- `n_datasets_ref` and `fraction_datasets_ref` for event ubiquity
- `splice_event_ubiquity_by_dataset_human_v1.tsv.gz` for runtime same-dataset exclusion
- per-gene `n_unique_event_groups_ref`
- per-gene `n_studies_ref`
- per-gene `n_studies_high_confidence_ref`
- per-gene `fraction_low_confidence_events_ref`
- per-gene `median_unique_groups_per_study`
- `splice_gene_event_burden_by_dataset_human_v1.tsv.gz` for runtime same-dataset exclusion

This keeps one very large source study from defining the prior on its own.

## Same-dataset exclusion at runtime

When the target run can identify its own `source_dataset`, runtime prior use defaults to:

- `--bundle_same_dataset_policy exclude`

Under `exclude`, matching source-dataset contributions are removed from:

- event ubiquity priors
- study-count recurrence used by conservative impact modulation
- gene-burden reference counts

If no non-self reference information remains after exclusion, the relevant prior becomes neutral and the run summary records that fallback explicitly.

Other policies:

- `warn`: keep the pooled prior but emit an explicit warning
- `fail`: refuse to run if matching source-dataset contributions are present
- `ignore`: use the pooled prior without warning

This is meant to protect against scoring `BRCA` with a bundle built from that same staged `BRCA` cohort by accident.

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
