# ATAC-seq to Gene Sets (Practical Guide)

`omics2geneset` ATAC converters map peak-level accessibility signals to
gene-level programs and gene sets. This guide is the practical ATAC entrypoint;
repository-wide framework guidance is in `README.md`.

## Install and Run

```bash
python -m pip install -U pip
python -m pip install -e ".[dev]"
omics2geneset list
```

Core CLI:

```bash
omics2geneset list
omics2geneset describe atac_bulk
omics2geneset describe atac_bulk_matrix
omics2geneset convert <converter_name> [flags]
omics2geneset validate <output_dir>
omics2geneset resources list
omics2geneset resources status --fast
omics2geneset resources status --check_schema --verify
omics2geneset resources describe ccre_ubiquity_hg19
omics2geneset resources describe ccre_ubiquity_hg38
omics2geneset resources fetch --preset atac_default_optional_hg19
omics2geneset resources fetch --preset atac_default_optional_hg38
omics2geneset resources manifest-validate
```

Quickstart (baseline bulk BED):

```bash
# Optional reference resources:
omics2geneset resources fetch --preset atac_default_optional_hg38
# omics2geneset resources fetch --preset atac_default_optional_hg19

omics2geneset convert atac_bulk \
  --peaks <peaks.bed.gz> \
  --gtf <genes.gtf.gz> \
  --organism human \
  --genome_build hg38 \
  --out_dir <out>

omics2geneset validate <out>
```

Resource catalog notes:

- Bundled manifest: `src/omics2geneset/resources/manifest.json`
- Default cache: `~/.cache/omics2geneset/resources` (override with `OMICS2GENESET_RESOURCES_DIR`)
- Preferred bundle workflow: extract a build-specific bundle (`hg19` or `hg38`), generate a local manifest with `--layout direct --genome-build <build>`, and pass `--resources_dir <bundle_root>` directly (no `resources fetch` required).
- `resources fetch` supports individual ids and presets.
- Human build presets: `atac_default_optional_hg19`, `atac_default_optional_hg38`.
- `resources fetch` skips resources with missing URL by default (status `manual`).
- `resources status --check_schema` validates known resource file schemas (beyond existence/checksum).
- Resource `url` fields can be HTTP(S), `file://`, or plain Unix filesystem paths.
- `--manifest_mode overlay` (default) merges a custom manifest with bundled entries; use `replace` to use only custom entries.
- Bundle setup instructions and split-bundle packaging: `docs/assays/atac/reference_bundle.md`.

## Output Contract

Program extractors write:

1. `geneset.tsv` (selected program)
2. `geneset.meta.json` (provenance, parameters, summary)

Optional:

3. `geneset.full.tsv` (full nonzero score table when `--emit_full true`)
4. `genesets.gmt` (one or more exported gene sets when `--emit_gmt true`)
5. `run_summary.json` and `run_summary.txt` (execution/QC summary)

GMT defaults favor cleaner symbols:

- `--gmt_require_symbol true` drops rows whose symbol is missing or Ensembl-like.
- `--gmt_biotype_allowlist protein_coding` keeps only protein-coding genes by default (when available).
- default ATAC preset is `--program_preset connectable` (alias: `default`) and emits only two recommended outputs:
  - `linked_activity` with `nearest_tss`
  - `distal_activity` with `distance_decay`
- default calibration policy is `--calibration_methods auto_prefer_ref_ubiquity_else_none`:
  - uses `ref_ubiquity_penalty` when resource data is available
  - otherwise falls back to `none`
- `--program_preset qc`, `--program_preset experimental`, and `--program_preset all` are opt-in for broader/non-recommended outputs.
- `--program_methods` controls only program families (link-level and family-level scoring), not contrast options.
- `--use_reference_bundle true` (default) allows reference-backed contrasts; `false` forces contrast behavior to `none`.
- `--use_reference_bundle false` is the explicit opt-out.
- default GMT selection plans are compact: `--gmt_topk_list 200` and `--gmt_mass_list ""`.
- default GMT size guardrails are strict: `--gmt_min_genes 100`, `--gmt_max_genes 500`, and `--emit_small_gene_sets false`.
- Human builds currently supported out of the box: `hg19` (`GRCh37`) and `hg38` (`GRCh38`).

Not recommended by default:

- `promoter_activity` and promoter-overlap linkage outputs for primary discovery (use mainly for QC)
- `enhancer_bias` and `tfidf_distal` as routine defaults (use in exploratory/experimental settings)
- `atlas_residual` (resource-backed and sensitive to score-definition matching)

`geneset.tsv` columns:

- required: `gene_id`, `score`, `rank`
- optional: `weight`, `gene_symbol`

For ATAC defaults (`--normalize within_set_l1`), `weight` is normalized only within selected genes.
Default ATAC usage attempts reference-backed methods when bundle resources are available.
If resources are missing, methods are skipped by default (`--resource_policy skip`), and converters print runtime warnings for each skipped contrast with hints on how to enable it.

### Why you might see fewer gene sets now

- `omics2geneset` now skips GMT entries smaller than `--gmt_min_genes` by default.
- For `condition_within_group`, groups are skipped when they fail cell and donor support gates (`--min_cells_per_group`, `--min_cells_per_condition`, `--min_donors_per_condition`) or dataset-level donor coverage (`--min_total_donors_per_condition`).
- Use `--emit_small_gene_sets true` only when you intentionally want small sets for debugging or toy runs.

### Implementation details: exact skip decision flow

For all ATAC converters (`atac_bulk`, `atac_bulk_matrix`, `atac_sc_10x`), GMT emission uses the same shared decision logic:

1. Filter candidate genes by optional `--gmt_biotype_allowlist`.
2. For each output variant (primary, or `pos/neg` when signed splitting is enabled), keep only rows with `score > 0`.
3. Build unique gene tokens using `--gmt_prefer_symbol` / `--gmt_require_symbol` (duplicates are dropped while preserving order).
4. If no positive-token genes remain, skip and warn.
5. If positive-token count is `< --gmt_min_genes`:
   - skip and warn when `--emit_small_gene_sets false` (default)
   - emit and warn when `--emit_small_gene_sets true`
6. If positive-token count is between `--gmt_min_genes` and `--gmt_min_genes + 50`, emit a marginal-signal warning.
7. For each plan in `--gmt_topk_list` / `--gmt_mass_list`, create a set from the positive-score rows only.
8. After token de-duplication, if the resulting set falls below `--gmt_min_genes`, apply the same skip/emit rule from step 5.

All decisions are written to metadata:

- `gmt.requested_outputs`
- `gmt.emitted_outputs`
- `gmt.skipped_outputs`
- `gmt.diagnostics`

## Conceptual model and mapping to atac-seq_methods.tex

Shared ATAC pipeline stages:

1. Define a peak statistic `x_p` from input peaks/matrices.
2. Transform to nonnegative peak weights `alpha_p = phi(x_p)`.
3. Apply optional external calibration (`--calibration_methods`) after `x_p` is defined.
4. Link peaks to genes with `L_pg` (`--link_method`) and compute gene scores `s_g = sum_p alpha_p * L_pg`.
5. Apply program-family views (`--program_preset` / `--program_methods`).
6. Select compact gene sets (`--select`, `--top_k`) and optional normalized weights (`--normalize`), then export GMT.

Important distinction:

- `--study_contrast` is the study-design contrast that defines `x_p` (for example, `group_vs_rest`, `condition_within_group`, or baseline/no contrast).
- `--calibration_methods` is external calibration applied after `x_p` is defined (`none`, `ref_ubiquity_penalty`, `atlas_residual`).
- Backward-compatible aliases remain available with deprecation warnings: `--contrast` -> `--study_contrast`, `--contrast_methods` -> `--calibration_methods`.

CLI to theory crosswalk:

| CLI flags (guide) | Model object (methods) | Defined in this note |
|---|---|---|
| `--peaks`, `--peak_matrix_tsv`, `--matrix_dir` | Peak coordinates and peak statistic `x_p` | Methods: Mental model (`sec:atac_mental_model`), Peak weights (`sec:peak_weights`) |
| `--gtf`, `--organism`, `--genome_build` | Gene coordinates (TSS/promoter/locus) | Methods: Linkage models (`sec:linkage_models`) |
| `--study_contrast` and related condition flags | Study-design contrast defining `x_p` | Methods: Study-design contrasts (`sec:study_design_contrast`) |
| `--peak_weight_transform` | Transform `alpha_p = phi(x_p)` | Methods: Peak weights (`sec:peak_weights`) |
| `--calibration_methods` | External calibration (peak-level or gene-level) | Methods: Calibration-method axis (`sec:contrast_axis`), `eq:ref_idf`, `eq:atlas_logratio`, `eq:atlas_z` |
| `--ref_ubiquity_resource_id`, `--atlas_resource_id`, `--atlas_metric`, `--atlas_min_raw_quantile`, `--atlas_use_log1p` | Reference resources and stabilization choices | Methods: Resources (`sec:resources`), M3.1/M3.2 |
| `--link_method` | Linkage model `L_pg` | Methods: Linkage models (`sec:linkage_models`) |
| `--program_preset`, `--program_methods` | Program-family view `Psi_m` | Methods: Program-family axis (`sec:program_axis`), catalog (`sec:program_catalog`) |
| `--select`, `--top_k`, `--normalize` | Gene set extraction operator and normalization | Methods: Set extraction (`sec:set_extraction`), `eq:within_program_l1` |
| `--emit_gmt`, `--gmt_topk_list`, `--gmt_min_genes`, `--gmt_max_genes` | GMT export of selected gene sets | Methods: GMT export (`sec:gmt_export`) |

## Extending omics2geneset

### Add a new converter in omics2geneset

1. Add `src/omics2geneset/converters/<name>.py` with `run(args) -> dict`.
2. Add CLI parser wiring in `src/omics2geneset/cli.py`.
3. Register converter in `src/omics2geneset/registry.py`.
4. Add converter spec JSON under `src/omics2geneset/converters/specs/<name>.json`.
5. Emit the contract files (`geneset.tsv`, `geneset.meta.json`, optional `geneset.full.tsv`).
6. Add toy fixtures and tests in `tests/`.

### Add a new mode to an existing omics2geneset converter

1. Add CLI flags.
2. Implement mode logic in converter/core modules.
3. Record resolved algorithmic parameters in metadata.
4. Add tests for defaults, edge cases, and validation behavior.

## Extractor: atac_bulk

ATAC program generation follows three independent axes:

- linkage axis (`--link_method`): how peaks are linked to genes
- calibration axis (`--calibration_methods`): how peak/gene signals are calibrated (`none`, `ref_ubiquity_penalty`, `atlas_residual`)
- program-family axis (`--program_methods` / `--program_preset`): how linked scores are turned into biologically motivated program views (`linked_activity`, `promoter_activity`, `distal_activity`, `enhancer_bias`, and `tfidf_distal` for scATAC)

### Theory cross-reference

- Methods: Peak weights (`sec:peak_weights`)
- Methods: Linkage models (`sec:linkage_models`)
- Methods: Calibration methods (`sec:contrast_axis`)
- Methods: Program-family axis and catalog (`sec:program_axis`, `sec:program_catalog`)
- Methods: Set extraction and GMT export (`sec:set_extraction`, `sec:gmt_export`)

Recommended default (`--program_preset connectable`) is intentionally not a full cross-product. It emits:

- `linked_activity` at `link_method=nearest_tss`
- `distal_activity` at `link_method=distance_decay`
- one calibration outcome per pair via `auto_prefer_ref_ubiquity_else_none`

Use `--program_preset all` for full method-development cross-product behavior.

### Required inputs

- `--peaks`: BED/narrowPeak-like intervals (`.gz` supported)
- `--gtf`: gene annotation
- `--out_dir`
- `--organism`
- `--genome_build`

Optional peak weights:

- `--peaks_weight_column` (default `5`; 1-based indexing where column 1 is `chrom`), or
- `--peak_weights_tsv` with `chrom,start,end,weight`

### Extraction modes (concept + scientific intent)

- Peak-to-gene linking (`--link_method`):
  - `promoter_overlap`: promoter-proximal assignment
  - `nearest_tss`: nearest-gene assignment
  - `distance_decay`: multi-gene weighted assignment by distance
  - `external`: user-provided region->gene linkage table (ABC/Hi-C-style)
  - modes can be combined (for GMT) as comma-separated tokens, e.g. `all,external`
- Peak transform (`--peak_weight_transform`): `signed`, `abs`, `positive`, `negative`
- Program families:
- `--program_preset {none,default,connectable,qc,experimental,all}` controls ATAC program families.
- `--program_methods` allows explicit family override (`linked_activity,promoter_activity,distal_activity,enhancer_bias`).
- in `connectable/default`, only recommended pairs are emitted (linked+nearest, distal+distance).
- use `all` for full cross-product.
- Calibration methods:
  - `--calibration_methods {auto_prefer_ref_ubiquity_else_none,all,none,ref_ubiquity_penalty,atlas_residual}`
  - default is `auto_prefer_ref_ubiquity_else_none`
  - evaluated independently across each selected `--link_method`
  - `--use_reference_bundle {true,false}` defaults to `true`; set `false` to force `calibration_methods=none`.
  - `--resource_policy {skip,fail}` controls behavior when resource-backed methods cannot load resources. Recommended: `fail` for production runs.
  - `--ref_ubiquity_resource_id` and `--atlas_resource_id` select catalog resources.
  - `--atlas_metric {logratio,zscore}` (default `zscore`) plus `--atlas_min_raw_quantile` (default `0.95`) and `--atlas_use_log1p` stabilize atlas residual scoring.
- Program selection (`--select`): `top_k`, `quantile`, `threshold`, or `none`
- Normalization (`--normalize`):
  - `within_set_l1` (default): normalize only selected genes
  - `l1`: legacy global normalization
  - `none`: raw selected scores as weights
- Optional marker QC:
  - `--qc_marker_genes_tsv <path>` loads marker symbols/IDs and records top-rank hit counts in metadata/run summary.

### Outputs

- `out_dir/geneset.tsv` (selected program)
- `out_dir/geneset.full.tsv` (optional full nonzero table)
- `out_dir/genesets.gmt` (optional GMT export; default on)
- `out_dir/geneset.meta.json`
- `out_dir/run_summary.json` and `out_dir/run_summary.txt`

### Quickstart

```bash
omics2geneset convert atac_bulk \
  --peaks tests/data/toy_peaks.bed \
  --gtf tests/data/toy.gtf \
  --out_dir tests/tmp/readme_bulk \
  --organism human \
  --genome_build hg38 \
  --peak_weights_tsv tests/data/toy_peak_weights.tsv \
  --program_preset connectable \
  --calibration_methods auto_prefer_ref_ubiquity_else_none \
  --select top_k \
  --top_k 200 \
  --normalize within_set_l1
```

### Common pitfalls

- Genome build mismatch between peaks, `--gtf`, and reference resources.
- Wrong `--peaks_weight_column` (for example selecting rank/length instead of accessibility signal).
- Missing resources causing `--calibration_methods` to fall back under `--resource_policy skip`.
- Promoter-heavy outputs are often generic for baseline single-sample data; prefer distal programs for discovery.
- Check `run_summary.txt` to confirm which methods ran versus skipped.

## Extractor: atac_sc_10x

### Required inputs

- `--matrix_dir` containing:
  - `matrix.mtx` or `matrix.mtx.gz`
  - `barcodes.tsv` or `barcodes.tsv.gz`
  - `peaks.bed` or `features.tsv(.gz)` with genomic peak coordinates
- `--gtf`
- `--out_dir`
- `--organism`
- `--genome_build`

Optional:

- `--groups_tsv` (`barcode`, `group`) for per-group programs
- `--cell_metadata_tsv` (`barcode` plus metadata columns) to enable condition contrasts

### Theory cross-reference

- Methods: Study-design contrasts (`sec:study_design_contrast`)
- Methods: scATAC peak summaries (`sec:sc_peak_summaries`)
- Methods: TF-IDF equation (`eq:tfidf`)
- Methods: Linkage/program/calibration axes (`sec:linkage_models`, `sec:contrast_axis`, `sec:program_axis`)
- Methods: Set extraction and GMT export (`sec:set_extraction`, `sec:gmt_export`)

### Extraction modes (concept + scientific intent)

- Peak summary (`--peak_summary`): `sum_counts`, `mean_counts`, `frac_cells_nonzero`
- Study contrast (`--study_contrast`):
  - `group_vs_rest` (default when groups are provided) to capture group-specific accessibility
  - `condition_within_group` for per-group condition A vs B contrasts with OPEN/CLOSE programs
  - `none` for non-differential summaries
- Contrast metric (`--contrast_metric`): `log2fc` or `diff`
- Pseudocount (`--contrast_pseudocount`):
  - auto defaults: `1e-3` for `frac_cells_nonzero`, `1.0` for count-based summaries
- Condition metadata flags:
  - `--condition_column <name>` in `cell_metadata_tsv`
  - `--donor_column <name>` in `cell_metadata_tsv` (default `donor`)
  - `--condition_a <label>` and `--condition_b <label>` (optional when exactly two levels are present)
  - `--min_cells_per_group <int>` minimum cells required for a group to emit outputs (default `100`)
  - `--min_cells_per_condition <int>` minimum cells per condition within each group
  - `--min_donors_per_condition <int>` minimum donors per condition within each group (default `3`)
  - `--min_total_donors_per_condition <int>` minimum donors per condition overall (default `3`)
  - `--cell_imbalance_warn_ratio <float>` warn on extreme case/control cell imbalance (default `5.0`)
  - `--baseline_carry_through_corr_warn <float>` warn when condition contrast tracks baseline too strongly (default `0.9`)
- Peak transform (`--peak_weight_transform`):
  - default `positive` for opening programs
  - use `negative` for closing programs
- Program families:
  - `--program_preset {none,default,connectable,qc,experimental,all}` controls ATAC program families.
- `--program_methods` overrides with explicit families; scATAC supports `tfidf_distal` in addition to bulk methods.
  - in `connectable/default`, only recommended pairs are emitted (linked+nearest, distal+distance).
  - use `all` for full cross-product.
- Calibration methods:
  - `--calibration_methods {auto_prefer_ref_ubiquity_else_none,all,none,ref_ubiquity_penalty,atlas_residual}`
  - default is `auto_prefer_ref_ubiquity_else_none`
  - for `--study_contrast condition_within_group`, default/auto policy is forced to `none`; to run reference calibration explicitly, pass `--calibration_methods ref_ubiquity_penalty` (or another explicit non-auto choice)
  - evaluated independently across each selected `--link_method`
  - `--use_reference_bundle {true,false}` defaults to `true`; set `false` to force `calibration_methods=none`.
  - `--resource_policy {skip,fail}` controls behavior when resource-backed methods cannot load resources. Recommended: `fail` for production runs.
  - `--ref_ubiquity_resource_id` and `--atlas_resource_id` select catalog resources.
  - `--atlas_metric {logratio,zscore}` (default `zscore`) plus `--atlas_min_raw_quantile` (default `0.95`) and `--atlas_use_log1p` stabilize atlas residual scoring.
- Gene program selection and normalization: same controls as `atac_bulk`
- GMT guardrails:
  - by default, sets with fewer than `--gmt_min_genes` are skipped with warnings
  - set `--emit_small_gene_sets true` to allow emitting small sets

### Group/contrast skip behavior (implementation)

For `atac_sc_10x`, group-level gating is applied before GMT extraction:

1. For `--study_contrast condition_within_group`, if `--calibration_methods` is auto/default, runtime forces `calibration_methods=none` unless you explicitly request a non-auto method.
2. For `--study_contrast condition_within_group`, if condition metadata is present:
   - print per-group case/control cell and donor counts
   - skip groups with `n_cells < --min_cells_per_group`
   - skip groups when either condition has `< --min_cells_per_condition` cells
   - skip groups when either condition has `< --min_donors_per_condition` donors
   - warn when case/control cell imbalance exceeds `--cell_imbalance_warn_ratio`
   - optionally warn on high baseline carry-through using `--baseline_carry_through_corr_warn`
3. For `--study_contrast condition_within_group`, if overall donor coverage is too low (`--min_total_donors_per_condition`), explicit runs fail; implicit preset-selected runs fall back.
4. If no groups remain after condition QC filtering, conversion exits with an error.
5. If `condition_within_group` is selected implicitly by preset (not explicitly set by `--study_contrast`) but condition metadata is missing or unusable, the converter warns and falls back to:
   - `group_vs_rest` when groups are available
   - `none` when groups are not available
6. If `--study_contrast condition_within_group` is explicitly requested and required condition metadata is missing/unusable, conversion errors instead of falling back.

External linkage TSV format:

- required columns: `chrom`, `start`, `end`, `gene_id`, `link_weight`
- 0-based half-open coordinates
- `gene_id` is resolved against GTF IDs; unversioned Ensembl IDs are accepted

Resource-backed table formats:

- `ref_ubiquity_penalty` resource table:
  - required: `chrom`, `start`, `end`, and either `idf_ref` or (`df_ref`, `n_ref`)
- `atlas_residual` resource table:
  - required: `gene_id`, `median_score`, `mad_score`
  - accepted aliases: `median`, `mad`

Common resource IDs:

- `ccre_ubiquity_hg19`, `ccre_ubiquity_hg38` for `ref_ubiquity_penalty` (human)
- `atac_reference_profiles_hg19`, `atac_reference_profiles_hg38` for `atlas_residual` (human)
- `ccre_ubiquity_mm10`, `atac_reference_profiles_mm10` remain optional legacy entries

### Outputs

Without groups:

- `out_dir/geneset.tsv`
- `out_dir/geneset.full.tsv` (optional)
- `out_dir/genesets.gmt` (optional; default on)
- `out_dir/geneset.meta.json`
- `out_dir/run_summary.json`
- `out_dir/run_summary.txt`

With groups:

- `out_dir/group=<GROUP>/geneset.tsv`
- `out_dir/group=<GROUP>/geneset.full.tsv` (optional)
- `out_dir/group=<GROUP>/genesets.gmt` (optional; default on)
- `out_dir/group=<GROUP>/geneset.meta.json`
- `out_dir/group=<GROUP>/run_summary.json`
- `out_dir/group=<GROUP>/run_summary.txt`
- `out_dir/genesets.gmt` (combined all-group GMT file)
- `out_dir/manifest.tsv`

Top-level grouped validation:

```bash
omics2geneset validate <out_dir>
```

### Quickstart (cluster-specific programs)

```bash
omics2geneset convert atac_sc_10x \
  --matrix_dir tests/data/toy_10x_mtx \
  --gtf tests/data/toy.gtf \
  --groups_tsv tests/data/barcode_groups.tsv \
  --out_dir tests/tmp/readme_sc \
  --organism human \
  --genome_build hg38 \
  --program_preset connectable \
  --study_contrast group_vs_rest \
  --select top_k \
  --top_k 200 \
  --normalize within_set_l1
```

### Quickstart (connectable condition-within-group OPEN/CLOSE)

```bash
omics2geneset convert atac_sc_10x \
  --matrix_dir tests/data/toy_10x_mtx \
  --gtf tests/data/toy.gtf \
  --groups_tsv tests/data/barcode_groups.tsv \
  --cell_metadata_tsv tests/data/barcode_conditions.tsv \
  --condition_column condition \
  --condition_a treated \
  --condition_b control \
  --min_cells_per_condition 1 \
  --out_dir tests/tmp/readme_sc_connectable \
  --organism human \
  --genome_build hg38 \
  --program_preset connectable
```

## Extractor: atac_bulk_matrix

Use this converter for condition contrasts across bulk ATAC samples using a peak-by-sample matrix.

### Theory cross-reference

- Methods: Study-design contrasts (`sec:study_design_contrast`)
- Methods: Peak weights and transforms (`sec:peak_weights`)
- Methods: Linkage and calibration axes (`sec:linkage_models`, `sec:contrast_axis`)
- Methods: Program catalog and defaults (`sec:program_catalog`, `sec:presets`)
- Methods: Set extraction and GMT export (`sec:set_extraction`, `sec:gmt_export`)

Required inputs:

- `--peak_matrix_tsv`: rows are peaks, columns are samples (numeric values)
- `--peak_bed`: row-aligned peak coordinates
- `--sample_metadata_tsv`: includes `sample_id` + `condition` columns (or overrides)
- `--gtf`, `--organism`, `--genome_build`, `--out_dir`

Key study-contrast and calibration flags:

- `--condition_a`, `--condition_b`
- `--contrast_metric {log2fc,diff}`
- `--contrast_pseudocount` (used for `log2fc`)
- `--calibration_methods {auto_prefer_ref_ubiquity_else_none,all,none,ref_ubiquity_penalty,atlas_residual}` (default `auto_prefer_ref_ubiquity_else_none`)
- `--use_reference_bundle {true,false}` defaults to `true`; set `false` to force reference-backed contrasts off.
- `--resource_policy {skip,fail}` recommended as `fail` for production runs.

This converter emits direction-aware outputs:

- OPEN programs: peaks more accessible in condition A
- CLOSE programs: peaks less accessible in condition A

`genesets.gmt` includes both directions with semantic names (assay, dataset label, condition labels, direction, program method).

Output directory also includes `run_summary.json` and `run_summary.txt` with link/contrast execution details and optional marker QC.

Quickstart:

```bash
omics2geneset convert atac_bulk_matrix \
  --peak_matrix_tsv tests/data/toy_bulk_peak_matrix.tsv \
  --peak_bed tests/data/toy_peaks.bed \
  --sample_metadata_tsv tests/data/toy_bulk_sample_metadata.tsv \
  --condition_a case \
  --condition_b control \
  --gtf tests/data/toy.gtf \
  --out_dir tests/tmp/readme_bulk_matrix \
  --organism human \
  --genome_build hg38 \
  --program_preset connectable \
  --gmt_min_genes 1 \
  --gmt_max_genes 10 \
  --gmt_topk_list 3 \
  --gmt_mass_list ""
```
