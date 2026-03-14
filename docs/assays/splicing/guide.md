# Alternative Splicing Guide

This assay family adds two converter entrypoints plus two workflow helpers:

- `splice_event_diff`: event-level differential splicing table -> signed gene program and GMT sets
- `splice_event_matrix`: PSI matrix plus sample metadata -> one or more contrast-specific gene programs
- `workflows splice_prepare_public`: narrow public-data normalizer for TCGA SpliceSeq-like PSI matrices
- `workflows splice_prepare_reference_bundle`: compact alias/ubiquity/impact/burden bundle builder

The modeling stance is deliberately narrow:

- event-to-gene mapping, not full transcript reconstruction
- uncertainty-aware event scoring
- conservative impact priors that shrink toward neutral
- aggregation that resists genes winning only because they have many measurable events

## Quickstart

### 1. Optional: build a local splicing bundle

```bash
geneset-extractors workflows splice_prepare_reference_bundle \
  --sources_tsv <sources.tsv> \
  --out_dir <splice_bundle_dir> \
  --organism human \
  --bundle_id splice_human_v1
```

This writes compact files such as:

- `splice_event_aliases_human_v1.tsv.gz`
- `splice_event_ubiquity_human_v1.tsv.gz`
- `splice_event_impact_human_v1.tsv.gz`
- `splice_gene_event_burden_human_v1.tsv.gz`
- `bundle_provenance.json`
- `local_resources_manifest.json`

If you point `--resources_dir` at that directory, both splicing converters auto-resolve the four splicing resources by default for human runs.

### 2. Convert an event-level differential splicing table

```bash
geneset-extractors convert splice_event_diff \
  --splice_tsv <splice_event_diff.tsv> \
  --out_dir <out_dir> \
  --organism human \
  --genome_build hg38 \
  --resources_dir <splice_bundle_dir>
```

Minimal no-bundle run:

```bash
geneset-extractors convert splice_event_diff \
  --splice_tsv <splice_event_diff.tsv> \
  --out_dir <out_dir> \
  --organism human \
  --genome_build hg38 \
  --use_reference_bundle false
```

Useful explicit flags:

- `--tool_family auto|generic|leafcutter|majiq|whippet|tcga_spliceseq`
- `--score_mode auto|stat|delta_psi_times_neglog10p|delta_psi_times_confidence|custom_column`
- `--confidence_weight_mode none|pvalue|probability|read_support|combined`
- `--delta_psi_soft_floor 0.05`
- `--delta_psi_soft_floor_mode auto|off`
- `--impact_mode none|conservative|custom_bundle`
- `--event_dup_policy highest_confidence|max_abs|mean|sum`
- `--gene_burden_penalty_mode none|current_input|reference_bundle|auto`
- `--ambiguous_gene_policy drop|split_equal|first`

LeafCutter nuance:

```bash
--tool_family leafcutter --cluster_stats_tsv <cluster_stats.tsv>
```

If row-level p-values are missing, the converter will merge cluster significance onto intron rows when possible.

### 3. Convert a PSI matrix with within-study contrasts

```bash
geneset-extractors convert splice_event_matrix \
  --psi_matrix_tsv <psi_matrix.tsv> \
  --sample_metadata_tsv <sample_metadata.tsv> \
  --event_metadata_tsv <event_metadata.tsv> \
  --study_contrast condition_a_vs_b \
  --condition_a case \
  --condition_b control \
  --out_dir <out_dir> \
  --organism human \
  --genome_build hg38 \
  --resources_dir <splice_bundle_dir>
```

Grouped contrasts are also supported:

- `--study_contrast group_vs_rest`
- `--study_contrast condition_within_group`
- `--study_contrast baseline`

The matrix front-end computes standardized event-level contrasts and then reuses the same scorer as `splice_event_diff`.

### 4. Stage TCGA SpliceSeq-like public PSI data

```bash
geneset-extractors workflows splice_prepare_public \
  --input_mode tcga_spliceseq \
  --psi_tsv <tcga_spliceseq.tsv> \
  --sample_annotations_tsv <sample_annotations.tsv> \
  --out_dir <prepared_dir> \
  --organism human \
  --genome_build hg38 \
  --study_id <study_id> \
  --study_label <study_label>
```

This writes:

- `psi_matrix.tsv`
- `sample_metadata.tsv`
- `event_metadata.tsv`
- `sample_id_map.tsv`
- `event_id_map.tsv`
- `prepare_summary.json`
- `bundle_source_row.tsv`

`bundle_source_row.tsv` is deliberately shaped so many staged studies can be concatenated and then referenced from `--sources_tsv` for `splice_prepare_reference_bundle`.

For raw MD Anderson TCGA SpliceSeq downloads, the staging workflow now accepts the common public header aliases directly:

- `symbol` as `gene_symbol`
- `as_id` as `event_id`

Preferred sample metadata path:

- pass an explicit `--sample_annotations_tsv`

Fallback path when that table is unavailable:

- columns ending in `_Norm` are treated as adjacent-normal
- unsuffixed sample columns are treated as tumor
- the standardized `sample_id` is the header with `_Norm` stripped
- if the standardized sample ids would collide, the workflow raises and asks for explicit sample annotations

`prepare_summary.json` records whether metadata were explicit or inferred, the exact inference rule, tumor versus adjacent-normal counts, and whether TCGA barcode validation was attempted.

### 5. Validate outputs

```bash
geneset-extractors validate <out_dir>
```

## Input contracts

### `splice_event_diff`

Expected table shape:

- one row per event or event summary
- signed effect columns such as `stat` or `delta_psi`
- optional significance / posterior / read-support columns
- event-to-gene mapping columns
- optional coordinate and annotation columns

Common logical columns:

- `event_id`
- `event_group`
- `event_type`
- `gene_id`
- `gene_symbol`
- `chrom`, `start`, `end`, `strand`
- `stat`
- `delta_psi`
- `padj` / `pvalue`
- `probability`
- `read_support`
- `novel_flag`
- `annotation_status`

### `splice_event_matrix`

Expected files:

- `psi_matrix.tsv`: one row per event, one column per sample
- `sample_metadata.tsv`: at least `sample_id`, and usually `condition` and/or `group`
- optional `event_metadata.tsv`
- optional `coverage_matrix.tsv`

## Defaults

Default direct-converter behavior:

- `score_mode=auto`
- `score_transform=signed`
- `confidence_weight_mode=combined`
- `impact_mode=conservative`
- `event_dup_policy=highest_confidence`
- `gene_aggregation=signed_topk_mean`
- `gene_topk_events=3`
- `gene_burden_penalty_mode=auto`
- `delta_psi_soft_floor=0.05`
- `delta_psi_soft_floor_mode=auto`
- `ambiguous_gene_policy=drop`
- `use_reference_bundle=true`
- `resource_policy=skip`
- `select=top_k`
- `top_k=200`
- `normalize=within_set_l1`
- `emit_gmt=true`
- `gmt_split_signed=true`

Default matrix behavior:

- `effect_metric=welch_t`
- `missing_value_policy=min_present`
- `min_samples_per_condition=3`
- `min_present_per_condition=3`

## Warnings and missing-bundle behavior

When `resource_policy=skip` and bundle resources are missing, the converters:

- warn to stderr
- record missing resources in metadata and run summaries
- continue with neutral priors

The run summaries also surface interpretable QC such as:

- event-type composition
- canonicalization confidence counts
- retained-intron fraction
- novel-event fraction
- low-support fraction
- single-event-gene fraction
- delta-PSI soft-floor downweight counts
- selected-event low-confidence prior matches

The main warning heuristics are:

- retained intron dominance
- novel event dominance
- low confidence dominance
- likely event multiplicity bias

## How the new safeguards work

### Canonicalization confidence

The public workflow and bundle now distinguish:

- `canonicalization_status=coordinate_canonical`, `canonicalization_confidence=high`
- `canonicalization_status=raw_id_fallback`, `canonicalization_confidence=low`

Low-confidence keys are cohort-local conveniences, not globally harmonized splice identities. Runtime priors can still use them, but their effect is shrunk much more strongly toward neutral.

### Gene burden penalty

After event-to-gene aggregation, the splicing runtime applies a conservative multiplicity penalty:

- if a reference burden table is available, `gene_burden_penalty_mode=auto` uses it
- otherwise it falls back to the current input after duplicate collapse

This is meant to stop genes with many measurable events from winning on opportunity alone, not to erase them.

### Delta-PSI soft floor

When a usable `delta_psi` exists, very small absolute PSI shifts are softly downweighted instead of hard-dropped. This matters most in large cohorts where tiny but precise shifts can otherwise dominate significance-only scoring.
