# Morphology Profiling to Gene Sets (Practical Guide)

`geneset-extractors` morphology extraction consumes precomputed morphology profiles and compares them to a reference perturbation library. It does not compute image features from raw images.

Primary converter:

- `morphology_profile_query`

Primary workflow:

- `workflows jump_prepare_reference_bundle`

Theory and equations: `docs/assays/morphology/methods.tex`.
Reference bundle instructions: `docs/assays/morphology/reference_bundle.md`.

## What the extractor does

For each query morphology profile or query group, the extractor:

1. aligns shared numeric features between query and reference,
2. optionally standardizes using bundle feature statistics,
3. computes query-reference similarity,
4. routes compound matches through compound targets and genetic matches directly to genes,
5. balances compound and genetic modalities,
6. emits one or more gene-weighted programs and GMT sets.

Current specificity defaults are retrieval-oriented:

- prefer same-modality neighbors when query modality is known,
- use cross-modality neighbors mainly when they reinforce the same local target signal,
- treat `--max_reference_neighbors` as an upper bound rather than a fixed count,
- penalize genes that recur across many reference perturbations.

Current retrieval modes:

- `--mode direct_target`
  - strict target-first pooling
  - best when distributed same-target support exists
  - intentionally strict; on weak held-out queries it may emit very small outputs rather than broadening to a possibly wrong family
- `--mode mechanism`
  - family/mechanism-aware reranking and expansion
  - best when morphology is class-correct but not exact-target-correct
  - reaches its intended behavior only when the bundle or explicit inputs include `target_annotations`
  - chooses expansion labels from the retained mechanism neighborhood, not from reranked core genes
  - uses the hierarchy `pathway_seed -> target_class -> mechanism_label -> target_family`
- `--mode hybrid`
  - writes both `geneset.core.tsv` and `geneset.expanded.tsv`
  - keeps `geneset.tsv` as the preferred variant for downstream compatibility
  - computes the strict core branch and the mechanism branch independently
  - without target annotations, hybrid is mostly a strict core plus a conservative mechanism fallback rather than full family-aware expansion
  - for `orf` and `crispr` queries, expansion confidence still favors same-modality family/mechanism support
  - for `compound` queries, expansion may legitimately rely on coherent genetic neighbors and does not require same-modality support

Default interpretation:

- `--polarity similar` asks which known perturbations induce similar morphology.
- `--polarity opposite` asks which perturbations induce anti-correlated morphology.
- `--polarity both` requests both as separate programs, but opposite-polarity outputs may be suppressed if specificity is too low.

## Inputs

Two modes are supported.

### Explicit-file mode

Required:

- `--query_profiles_tsv`
- `--reference_profiles_tsv` or `--reference_profiles_parquet`
- `--reference_metadata_tsv`
- `--compound_targets_tsv`
- `--out_dir`
- `--organism`
- `--genome_build`

Optional:

- `--query_metadata_tsv`
- `--group_query_by`
- `--feature_schema_tsv`
- `--feature_stats_tsv`
- `--target_annotations_tsv`
  - overrides the packaged canonical annotation table used by the bundle builder

### Bundle-driven mode

Required:

- `--query_profiles_tsv`
- `--reference_bundle_id`
- `--resources_dir`
- `--out_dir`
- `--organism`
- `--genome_build`

Recommended bundle ids:

- `morphology_jump_target_pilot_u2os_48h_v1`
- `morphology_jump_target_pilot_a549_48h_v1`

## Quickstart: explicit files

```bash
geneset-extractors convert morphology_profile_query \
  --query_profiles_tsv path/to/query_profiles.tsv \
  --query_metadata_tsv path/to/query_metadata.tsv \
  --group_query_by query_group \
  --reference_profiles_tsv path/to/reference_profiles.tsv.gz \
  --reference_metadata_tsv path/to/reference_metadata.tsv.gz \
  --compound_targets_tsv path/to/compound_targets.tsv.gz \
  --feature_schema_tsv path/to/feature_schema.tsv.gz \
  --feature_stats_tsv path/to/feature_stats.tsv.gz \
  --out_dir results/morphology_query \
  --organism human \
  --genome_build hg38

geneset-extractors validate results/morphology_query
```

## Quickstart: bundle-driven

```bash
geneset-extractors convert morphology_profile_query \
  --query_profiles_tsv path/to/query_profiles.tsv \
  --query_metadata_tsv path/to/query_metadata.tsv \
  --group_query_by query_group \
  --reference_bundle_id morphology_jump_target_pilot_u2os_48h_v1 \
  --resources_dir /path/to/morphology_bundle_root \
  --out_dir results/morphology_query \
  --organism human \
  --genome_build hg38
```

If the bundle is missing and `--resource_policy skip`, the converter warns and stops only if required reference files cannot be resolved. If `--resource_policy fail`, missing bundles are fatal.

## Defaults

Default converter behavior is conservative and designed to produce interpretable outputs:

- `--mode direct_target`
- `--query_aggregate median`
- `--similarity_metric cosine`
- `--similarity_power 1.0`
- `--polarity similar`
- `--max_reference_neighbors 20`
- `--same_modality_first true`
- `--cross_modality_penalty 0.35`
- `--adaptive_neighbors true`
- `--mutual_neighbor_filter true`
- `--min_similarity 0.0`
- `--control_calibration mean_center`
- `--control_residual_components 2`
- `--control_min_profiles_for_residualization 5`
- `--hubness_penalty inverse_rank`
- `--gene_recurrence_penalty idf`
- `--min_specificity_confidence_to_emit_opposite medium`
- `--compound_weight 0.5 --genetic_weight 0.5`
- `--select top_k --top_k 200`
- `--normalize within_set_l1`
- `--emit_gmt true`
- `--gmt_topk_list 200`
- `--gmt_min_genes 100 --gmt_max_genes 500`

## Query grouping

If your query file has replicate wells or multiple profiles for the same biological unit, provide a metadata table and use:

```bash
--query_metadata_tsv path/to/query_metadata.tsv --group_query_by query_group
```

This aggregates queries before similarity scoring. Default aggregation is median.

## Building a reference bundle from JUMP-like data

The workflow command creates a local reusable bundle from well-level morphology profiles plus metadata:

```bash
geneset-extractors workflows jump_prepare_reference_bundle \
  --profile_paths path/to/profiles_plate1.tsv,path/to/profiles_plate2.tsv \
  --experimental_metadata_tsv path/to/experimental_metadata.tsv \
  --compound_targets_tsv path/to/jump_target_compound_metadata.tsv \
  --cell_type_filter U2OS \
  --timepoint_filter 48 \
  --require_same_timepoint_across_modalities true \
  --out_dir results/jump_u2os_48h_bundle
```

By default, the bundle builder now injects the packaged canonical morphology target-annotation table:

- `src/geneset_extractors/resources/morphology_target_annotations_human_v1.tsv`

Use:

```bash
--target_annotations_tsv path/to/custom_target_annotations.tsv
```

to override it, or:

```bash
--use_default_target_annotations false
```

to build an annotation-free bundle for strict fallback comparisons.

Minimal working bundle-to-converter path:

```bash
geneset-extractors workflows jump_prepare_reference_bundle \
  --profile_paths path/to/profiles.tsv \
  --experimental_metadata_tsv path/to/experimental_metadata.tsv \
  --compound_targets_tsv path/to/jump_target_compound_metadata.tsv \
  --cell_type_filter U2OS \
  --timepoint_filter 48 \
  --out_dir results/jump_u2os_48h_bundle

geneset-extractors convert morphology_profile_query \
  --query_profiles_tsv path/to/query_profiles.tsv \
  --query_metadata_tsv path/to/query_metadata.tsv \
  --group_query_by query_group \
  --reference_bundle_id morphology_jump_target_pilot_u2os_48h_v1 \
  --resources_dir results/jump_u2os_48h_bundle \
  --out_dir results/morphology_query \
  --organism human \
  --genome_build hg38
```

It writes:

- `reference_profiles.tsv.gz`
- `reference_metadata.tsv.gz`
- `compound_targets.tsv.gz`
- optional `target_annotations.tsv.gz`
- `feature_schema.tsv.gz`
- `feature_stats.tsv.gz`
- `<bundle_id>.bundle.json`

## Output layout

Grouped output layout mirrors other multi-program extractors:

- `manifest.tsv`
- `program=<QUERY>__polarity=<POLARITY>/geneset.tsv`
- optional `program=<QUERY>__polarity=<POLARITY>/geneset.core.tsv`
- optional `program=<QUERY>__polarity=<POLARITY>/geneset.expanded.tsv`
- `program=<QUERY>__polarity=<POLARITY>/geneset.meta.json`
- optional `program=<...>/geneset.full.tsv`
- optional `program=<...>/genesets.gmt`
- root `genesets.gmt` (combined)
- `run_summary.json` and `run_summary.txt`

## Routing diagnostics

When a morphology result looks surprising, inspect the per-program summary fields in:

- `program=<...>/run_summary.json`
- `program=<...>/geneset.meta.json`

Most useful fields:

- `raw_candidate_neighbor_ids` / `raw_candidate_neighbor_similarities`
  - top neighbors before hubness/context penalties and truncation
- `raw_candidate_neighbors_detail`
  - raw neighbor modality and top routed targets
- `retained_neighbors_detail`
  - neighbors that actually contributed after penalties/filtering
- `family_vote_summary_raw`
  - family-level vote totals in the pooled raw neighborhood
- `family_vote_summary_retained`
  - family-level vote totals after penalization/filtering
- `mechanism_vote_summary_raw`
- `mechanism_vote_summary_retained`
- `label_scores_raw`
- `label_scores_retained`
  - level-by-level support records for `pathway_seed`, `target_class`, `mechanism_label`, and `target_family`
- `expansion_decision`
  - machine-readable reason for whether expansion was allowed
  - includes `chosen_level`, `chosen_label`, `expansion_confidence`, `candidate_scope`, and whether a raw-vs-retained mismatch was observed
- `core_branch_neighbor_ids` / `mechanism_branch_neighbor_ids`
  - the strict exact-target branch and the broader mechanism branch are now reported separately
- `top_target_candidates`
  - routed target support before final gene-set extraction

Interpretation:

- if raw and retained label summaries disagree, that now lowers `expansion_confidence` instead of acting as a hard veto
- if `expansion_decision.reason` is `retained_label_support_too_weak`, the retained mechanism neighborhood never concentrated enough to justify expansion
- if `expansion_decision.chosen_level` is `target_class` or `pathway_seed`, the workflow found a narrower stable label and preferred it over the broader family
- for compound queries, expansion can still be valid when the strongest support comes from coherent genetic neighbors

## Common warnings and how to read them

- Low shared feature fraction:
  - query/reference feature spaces do not match well.
  - fix upstream feature selection or use a matched bundle.
- Reference controls excluded:
  - expected when bundle metadata marks control profiles.
- Control-derived calibration:
  - controls are not retrieval candidates, but they can still be used to estimate nuisance morphology structure before similarity scoring.
  - `--control_calibration mean_center` subtracts only the control mean.
  - `--control_calibration residualize_controls` subtracts the control mean and then projects query/reference vectors off the top control-derived nuisance axes.
  - if too few controls are available for residualization, the workflow falls back to mean-centering and records that fallback in summaries and metadata.
- Direct target versus mechanism mode:
  - `direct_target` pools positive evidence by target before hard neighbor truncation.
  - `mechanism` uses the full coherent retained neighborhood to choose the most specific supported annotation level, then expands locally from that mechanism branch.
  - `hybrid` emits both strict and expanded outputs and records which one is preferred; the strict branch does not define the expansion branch.
  - current expansion is confidence-weighted rather than hard-vetoed.
  - same-modality support still matters more for ORF/CRISPR queries than for compound queries.
  - bundle builds now include a canonical curated target-annotation table by default, so public/distributed bundles should normally be mechanism-ready.
  - some recurrent generic genes are intentionally left blank in that table to avoid noisy family expansion.
  - inspect `control_calibration`, branch-specific neighbor summaries, and `expansion_decision` if a result looks surprising.
- Many negative similarities ignored:
  - seen when `--polarity similar` but many anti-correlated matches exist.
  - rerun with `--polarity both` if you want both directions.
- Low retrieval confidence:
  - the retained morphology neighbors are weak or internally inconsistent.
  - check `run_summary.txt`, use a matched same-timepoint bundle, and inspect `top_neighbor_ids`.
- High retrieval confidence but low specificity confidence:
  - morphology neighbors are geometrically coherent, but the routed gene evidence is diffuse.
  - inspect `specificity_confidence`, `top10_gene_mass`, `neighbor_target_concentration`, `neighbor_primary_target_agreement`, `neighbor_top3_target_agreement`, and `high_hub_mass_fraction` in `geneset.meta.json`.
  - for opposite-polarity runs, low-specificity programs are suppressed by default unless you lower `--min_specificity_confidence_to_emit_opposite`.
- Generic recurrent genes dominate:
  - recurrent morphology genes such as receptor-family or stress-response genes can still appear.
  - default `--gene_recurrence_penalty idf` reduces this, but if you disable it you should expect broader outputs.
- Small or skipped GMT sets:
  - fewer than `--gmt_min_genes` genes survived selection.
  - for toy runs, use smaller `--gmt_min_genes` or `--emit_small_gene_sets true`.
  - warning format includes the actual gene count and threshold, for example:
    - `warning: small_gene_set_skipped program=Q1 polarity=similar n_genes=42 min_required=100`
- Artifact-family dominance:
  - top genes are dominated by broad receptor-family patterns.
  - advisory only; interpret cautiously and inspect the underlying matches.

## Recommended first run

For a first pass on Cell Painting data:

- use a bundle matched to the same cell line and timepoint,
- keep same-timepoint coherence unless you explicitly need a mixed-timepoint reference,
- aggregate replicate query wells,
- keep `--polarity similar`,
- treat `--polarity opposite` as experimental / secondary,
- keep the default hubness penalty on unless you are explicitly studying broad/generic morphological states,
- keep same-modality-first retrieval on when query modality metadata is available,
- keep the neighborhood narrow by default; increasing `--max_reference_neighbors` usually makes outputs broader and less specific,
- if controls capture a strong vehicle or plate axis, prefer `--control_calibration residualize_controls` over plain mean-centering,
- if you explicitly want exploratory opposite-polarity outputs, use `--min_specificity_confidence_to_emit_opposite low`,
- inspect `run_summary.txt` before over-interpreting GMT enrichments.

## Common mismatch: public JUMP contexts

Public JUMP pilot contexts do not always contain compound, ORF, and CRISPR perturbations at the same cell line and timepoint.

Default workflow behavior:

- do not silently mix timepoints,
- build a coherent same-timepoint bundle from the modalities that are actually available,
- warn explicitly about missing modalities.

If you intentionally want a mixed-timepoint bundle, opt in with:

```bash
--allow_mixed_timepoints true --require_same_timepoint_across_modalities false
```
