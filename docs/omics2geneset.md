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
```

## Output Contract

Program extractors write:

1. `geneset.tsv` (selected program)
2. `geneset.meta.json` (provenance, parameters, summary)

Optional:

3. `geneset.full.tsv` (full nonzero score table when `--emit_full true`)
4. `genesets.gmt` (one or more exported gene sets when `--emit_gmt true`)

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
- Peak transform (`--peak_weight_transform`): `signed`, `abs`, `positive`, `negative`
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
- Gene program selection and normalization: same controls as `atac_bulk`

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
