# omics2geneset Package Guide

`omics2geneset` is the first extractor family implemented in this repository.  
It converts assay-specific input signals into a common gene-program output contract.

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
omics2geneset convert <converter_name> [flags]
omics2geneset validate <output_dir>
omics2geneset resources list
omics2geneset resources status
omics2geneset resources fetch --preset atac_default_optional
```

Resource catalog notes:

- Bundled manifest: `src/omics2geneset/resources/manifest.json`
- Default cache: `~/.cache/omics2geneset/resources` (override with `OMICS2GENESET_RESOURCES_DIR`)
- `resources fetch` supports individual ids and presets.

## Output Contract

Program extractors write:

1. `geneset.tsv` (selected program)
2. `geneset.meta.json` (provenance, parameters, summary)

Optional:

3. `geneset.full.tsv` (full nonzero score table when `--emit_full true`)
4. `genesets.gmt` (one or more exported gene sets when `--emit_gmt true`)

GMT defaults favor cleaner symbols:

- `--gmt_require_symbol true` drops rows whose symbol is missing or Ensembl-like.
- `--gmt_biotype_allowlist protein_coding` keeps only protein-coding genes by default (when available).
- `--link_method all` runs promoter/nearest/distance linkage models for GMT generation by default.
- `--program_preset default` emits additional ATAC program families (promoter/distal/enhancer-bias; scATAC adds tfidf_distal).
- `--program_preset all` also attempts resource-backed methods (`ref_ubiquity_penalty`, `atlas_residual`).
- `--gmt_topk_list 100,200,500` and `--gmt_mass_list 0.5,0.8,0.9` emit six GMT sets per linkage model.

`geneset.tsv` columns:

- required: `gene_id`, `score`, `rank`
- optional: `weight`, `gene_symbol`

For ATAC defaults (`--normalize within_set_l1`), `weight` is normalized only within selected genes.

## Adding New Extractors or Modes

### Add a new extraction type (new converter)

1. Add `src/omics2geneset/converters/<name>.py` with `run(args) -> dict`.
2. Add CLI parser wiring in `src/omics2geneset/cli.py`.
3. Register converter in `src/omics2geneset/registry.py`.
4. Add converter spec JSON under `src/omics2geneset/converters/specs/<name>.json`.
5. Emit the contract files (`geneset.tsv`, `geneset.meta.json`, optional `geneset.full.tsv`).
6. Add toy fixtures and tests in `tests/`.

### Add a new extraction mode (existing converter)

1. Add CLI flags.
2. Implement mode logic in converter/core modules.
3. Record resolved algorithmic parameters in metadata.
4. Add tests for defaults, edge cases, and validation behavior.

## Extractor: atac_bulk

### Required inputs

- `--peaks`: BED/narrowPeak-like intervals (`.gz` supported)
- `--gtf`: gene annotation
- `--out_dir`
- `--organism`
- `--genome_build`

Optional peak weights:

- `--peaks_weight_column` (default `5`), or
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
  - `--program_preset {none,default,all}` controls additional ATAC program families for GMT output.
  - `--program_methods` allows explicit override (`linked_activity,promoter_activity,distal_activity,enhancer_bias,ref_ubiquity_penalty,atlas_residual`).
  - `--resource_policy {skip,fail}` controls behavior when resource-backed methods cannot load resources.
  - `--ref_ubiquity_resource_id` and `--atlas_resource_id` select catalog resources.
  - `--atlas_metric {logratio,zscore}` controls atlas residual scoring.
- Program selection (`--select`): `top_k`, `quantile`, `threshold`, or `none`
- Normalization (`--normalize`):
  - `within_set_l1` (default): normalize only selected genes
  - `l1`: legacy global normalization
  - `none`: raw selected scores as weights

### Outputs

- `out_dir/geneset.tsv` (selected program)
- `out_dir/geneset.full.tsv` (optional full nonzero table)
- `out_dir/genesets.gmt` (optional GMT export; default on)
- `out_dir/geneset.meta.json`

### Quickstart

```bash
omics2geneset convert atac_bulk \
  --peaks tests/data/toy_peaks.bed \
  --gtf tests/data/toy.gtf \
  --out_dir tests/tmp/readme_bulk \
  --organism human \
  --genome_build hg38 \
  --peak_weights_tsv tests/data/toy_peak_weights.tsv \
  --link_method promoter_overlap \
  --peak_weight_transform positive \
  --select top_k \
  --top_k 200 \
  --normalize within_set_l1
```

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

### Extraction modes (concept + scientific intent)

- Peak summary (`--peak_summary`): `sum_counts`, `mean_counts`, `frac_cells_nonzero`
- Contrast (`--contrast`):
  - `group_vs_rest` (default when groups are provided) to capture group-specific accessibility
  - `none` for non-differential summaries
- Contrast metric (`--contrast_metric`): `log2fc` or `diff`
- Pseudocount (`--contrast_pseudocount`):
  - auto defaults: `1e-3` for `frac_cells_nonzero`, `1.0` for count-based summaries
- Peak transform (`--peak_weight_transform`):
  - default `positive` for opening programs
  - use `negative` for closing programs
- Program families:
  - `--program_preset {none,default,all}` controls additional ATAC program families.
  - `--program_methods` overrides with explicit methods; scATAC supports `tfidf_distal` in addition to bulk methods.
  - `--resource_policy {skip,fail}` controls behavior when resource-backed methods cannot load resources.
  - `--ref_ubiquity_resource_id` and `--atlas_resource_id` select catalog resources.
  - `--atlas_metric {logratio,zscore}` controls atlas residual scoring.
- Gene program selection and normalization: same controls as `atac_bulk`

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

### Outputs

Without groups:

- `out_dir/geneset.tsv`
- `out_dir/geneset.full.tsv` (optional)
- `out_dir/genesets.gmt` (optional; default on)
- `out_dir/geneset.meta.json`

With groups:

- `out_dir/group=<GROUP>/geneset.tsv`
- `out_dir/group=<GROUP>/geneset.full.tsv` (optional)
- `out_dir/group=<GROUP>/genesets.gmt` (optional; default on)
- `out_dir/group=<GROUP>/geneset.meta.json`
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
  --contrast group_vs_rest \
  --contrast_metric log2fc \
  --peak_summary frac_cells_nonzero \
  --link_method promoter_overlap \
  --peak_weight_transform positive \
  --select top_k \
  --top_k 200 \
  --normalize within_set_l1
```
