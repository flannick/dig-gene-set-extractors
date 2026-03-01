# DNA Methylation to Gene Sets (Practical Guide)

`omics2geneset` methylation converters start from differential methylation tables (not raw IDATs) and emit:

- `geneset.tsv` (selected signed program)
- `geneset.meta.json` (provenance + parameters)
- optional `geneset.full.tsv`
- optional `genesets.gmt` (default directional `pos/neg` sets)
- `run_summary.json` and `run_summary.txt`

Theory and equations: `docs/methylation_methods.tex`.

## Converters

- `methylation_cpg_diff`: differential CpG/probe table to gene programs.
- `methylation_dmr_regions`: differential region table (DMRs) to gene programs.
- `methylation_dmr`: legacy converter for pre-aggregated gene-level DMR summaries.

## Quickstart: CpG/probe differential table

```bash
omics2geneset convert methylation_cpg_diff \
  --cpg_tsv path/to/cpg_diff.tsv \
  --probe_manifest_tsv path/to/probe_manifest.tsv \
  --gtf path/to/genes.gtf.gz \
  --organism human \
  --genome_build hg38 \
  --out_dir results/methylation_cpg

omics2geneset validate results/methylation_cpg
```

Required table columns:

- score inputs: `delta_beta` (or override `--delta_column`) and p-value/FDR column (`pvalue`/`padj`)
- coordinates, either:
  - `chrom` + `pos` (or `start`/`end`), or
  - `probe_id` with `--probe_manifest_tsv`

Probe manifest format:

- required: `probe_id`, `chrom`, and either `pos` or `start/end`
- default assumption: `pos` is 1-based (`--manifest_pos_is_0based false`)

Probe-only workflow with resources bundle defaults:

```bash
export OMICS2GENESET_RESOURCES_DIR=/path/to/resources_bundle

omics2geneset convert methylation_cpg_diff \
  --cpg_tsv path/to/cpg_probe_only.tsv \
  --gtf path/to/genes.gtf.gz \
  --organism human \
  --genome_build hg38 \
  --array_type 450k \
  --out_dir results/methylation_cpg
```

If `--probe_manifest_tsv` is omitted, the converter auto-tries a default manifest resource id from (`array_type`, `genome_build`) when resources are configured (`--resources_dir`, `--resources_manifest`, or `OMICS2GENESET_RESOURCES_DIR`).

## Quickstart: DMR region table

```bash
omics2geneset convert methylation_dmr_regions \
  --dmr_tsv path/to/dmr.tsv \
  --gtf path/to/genes.gtf.gz \
  --organism human \
  --genome_build hg38 \
  --out_dir results/methylation_dmr

omics2geneset validate results/methylation_dmr
```

Required DMR columns:

- `chrom`, `start`, `end`
- `delta_methylation` (or override `--delta_column`)
- `pvalue`/`padj` (for `score_mode=delta_times_neglog10p`)

## Default behavior

Defaults are conservative and directional:

- `--program_preset connectable` (default alias: `default`)
- core methods: `promoter_activity` and `distal_activity`
- preferred link pairings under `connectable`:
  - `promoter_activity` with `promoter_overlap`
  - `distal_activity` with `distance_decay`
- `--score_mode delta_times_neglog10p`
- activity-oriented signed scoring by default: contribution is proportional to `(-delta_methylation)`
- `--delta_orientation activity_oriented|raw` controls sign convention
- `--score_transform signed`
- `--aggregation weighted_mean` (reduces CpG/region density bias)
- `--select top_k --top_k 200`
- `--normalize within_set_l1`
- `--emit_gmt true --gmt_split_signed true`
- strict GMT guardrails: `--gmt_min_genes 100 --gmt_max_genes 500 --emit_small_gene_sets false`

## Link and program axes

- Link axis (`--link_method`): `promoter_overlap`, `nearest_tss`, `distance_decay`, `external`
- Program axis (`--program_preset` / `--program_methods`):
  - `promoter_activity`
  - `distal_activity`
  - `linked_activity` (mostly QC/secondary)

## Optional resource-backed inputs

### Probe manifest via resources

If you do not pass `--probe_manifest_tsv`, you can provide:

- `--probe_manifest_resource_id`
- `--resources_dir` (and optional `--resources_manifest`)

Placeholder resource IDs are included in `src/omics2geneset/resources/manifest.json`:

- `methylation_probe_manifest_450k_hg19`
- `methylation_probe_manifest_epic_hg19`
- `methylation_probe_manifest_450k_hg38`
- `methylation_probe_manifest_epic_hg38`

Runtime metadata and summaries record whether probe mapping came from an explicit TSV or a resolved resource id/path.

### Enhancer restriction for distal programs

Distal programs are non-promoter by default. To restrict to enhancer-overlapping regions:

- `--distal_mode enhancer_only`
- provide `--enhancer_bed` or `--enhancer_resource_id`

## Skip/fail policy for missing resources

- `--resource_policy skip` (default): warn and continue with methods that can run.
- `--resource_policy fail`: raise on missing required resources for requested methods.

## Common pitfalls

- Genome-build mismatch between methylation coordinates, GTF, and any manifests/resources.
- Probe tables with `probe_id` only but no manifest mapping.
- Wrong coordinate convention for manifest/input `pos` (1-based vs 0-based).
- Directionality interpretation:
  - with `--delta_orientation activity_oriented` (default):
    - `__pos` sets: hypomethylation-associated activity increase (`-delta > 0`)
    - `__neg` sets: hypermethylation-associated activity decrease (`-delta < 0`)
  - with `--delta_orientation raw`:
    - `__pos` sets: positive raw delta signal
    - `__neg` sets: negative raw delta signal
- The converter warns when:
  - sex chromosomes were dropped (`--drop_sex_chrom true`)
  - many rows fail chromosome/genome-build mapping
  - requested resources are missing under `--resource_policy skip`
  - emitted GMT sets show possible OR/KRT family dominance

## Optional artifact filtering

By default the tool only warns for possible family-dominant sets. To actively filter before selection/GMT:

- `--exclude_gene_symbol_regex` (repeatable; example: `--exclude_gene_symbol_regex '^OR'`)
- `--exclude_gene_symbols_tsv` (one symbol per line)

Filters run before top-k selection and GMT emission, so strong filters can reduce set size.

## Small DMR / toy-input guidance

Defaults are strict (`--gmt_min_genes 100`, `--emit_small_gene_sets false`), so very small DMR/CpG runs may emit no GMT sets.

Use one of:

- `--emit_small_gene_sets true` for debugging/toy runs
- lower `--gmt_min_genes` when biologically appropriate

## Documentation map

- package index: `docs/omics2geneset.md`
- ATAC guide: `docs/atac-seq2geneset.md`
- RNA guide: `docs/rna-seq2geneset.md`
- methylation methods note: `docs/methylation_methods.tex`
