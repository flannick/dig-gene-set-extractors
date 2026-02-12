# omics2geneset Package Guide

`omics2geneset` is the first implemented extractor family in this repository. It converts assay-specific upstream signals into standardized gene-weight outputs.

## What It Produces

All converters in this package write:

1. `geneset.tsv`
2. `geneset.meta.json`

`geneset.tsv` minimum columns:

- `gene_id`
- `weight`

Optional columns:

- `gene_symbol`
- `rank`

`geneset.meta.json` records schema version, converter parameters, input files + hashes, and summary metrics.

## Quick Start

From repo root:

```bash
python -m pip install -U pip
python -m pip install -e ".[dev]"
omics2geneset list
pytest -q
```

Core CLI commands:

```bash
omics2geneset list
omics2geneset describe atac_bulk
omics2geneset convert <converter_name> [flags]
omics2geneset validate <output_dir>
```

## Offline / Air-Gapped Install

For offline installation workflows (wheelhouse checks, bootstrap scripts, and no-build-isolation guidance), see:

- `docs/air_gapped_install.md`

If you only need to run package code directly from source:

```bash
PYTHONPATH=src python -m omics2geneset.cli list
```

## Package Layout

- `src/omics2geneset/cli.py`: CLI entrypoint
- `src/omics2geneset/registry.py`: converter registry + spec lookup
- `src/omics2geneset/converters/`: converter implementations
- `src/omics2geneset/converters/specs/`: machine-readable converter specs
- `src/omics2geneset/core/`: shared linking/scoring/normalization/metadata/validation
- `src/omics2geneset/io/`: assay input parsers
- `tests/`: unit tests + toy fixtures

## Adding New Extraction Types or Modes

### Add a new extraction type (new converter)

1. Implement a converter module in `src/omics2geneset/converters/` with `run(args)`.
2. Register it in `src/omics2geneset/registry.py`.
3. Add CLI wiring in `src/omics2geneset/cli.py`.
4. Add a converter spec JSON in `src/omics2geneset/converters/specs/<name>.json`.
5. Emit standard outputs (`geneset.tsv`, `geneset.meta.json`).
6. Add tests + toy data under `tests/`.

### Add a new extraction mode (within an existing converter)

1. Add/extend CLI flags.
2. Implement mode logic in converter/shared core.
3. Record resolved algorithmic parameters in metadata.
4. Add tests for defaults and mode branches.

## Extractor: `atac_bulk`

### Purpose

Convert bulk ATAC-seq peaks plus peak-level signal into gene-level weights.

### Required Inputs

- `--peaks PATH` (BED-like: `chrom`, `start`, `end`; gzipped supported)
- `--gtf PATH`
- `--out_dir PATH`
- `--organism {human,mouse}`
- `--genome_build STRING`

Peak weights can come from:

- `--peaks_weight_column INT` (default `5`), or
- `--peak_weights_tsv PATH` with `chrom,start,end,weight`

### Extraction Modes (Concept + Scientific Rationale)

Peak-to-gene linking (`--link_method`):

- `promoter_overlap`: promoter-proximal assignment
- `nearest_tss`: nearest-gene heuristic
- `distance_decay`: weighted multi-gene assignment by genomic distance

Peak transform (`--peak_weight_transform`):

- `signed`, `abs`, `positive`, `negative`

Normalization (`--normalize`):

- `l1` (default), `none`

### Outputs

- `out_dir/geneset.tsv`
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
  --peak_weight_transform abs \
  --normalize l1
```

## Extractor: `atac_sc_10x`

### Purpose

Convert 10x scATAC peak-by-cell matrices into gene-weight vectors (single output or grouped outputs).

### Required Inputs

- `--matrix_dir PATH` containing:
  - `matrix.mtx` or `matrix.mtx.gz`
  - `barcodes.tsv` or `barcodes.tsv.gz`
  - `peaks.bed` or `features.tsv(.gz)` coordinate rows
- `--gtf PATH`
- `--out_dir PATH`
- `--organism {human,mouse}`
- `--genome_build STRING`

Optional:

- `--groups_tsv PATH` (`barcode`, `group`)

### Extraction Modes (Concept + Scientific Rationale)

Group summarization (`--peak_summary`):

- `sum_counts`
- `mean_counts`
- `frac_cells_nonzero`

Then peak-to-gene linking (`promoter_overlap`, `nearest_tss`, `distance_decay`) and final normalization.

### Outputs

Without groups:

- `out_dir/geneset.tsv`
- `out_dir/geneset.meta.json`

With groups:

- `out_dir/group=<GROUP>/geneset.tsv`
- `out_dir/group=<GROUP>/geneset.meta.json`
- `out_dir/manifest.tsv`

### Quickstart

```bash
omics2geneset convert atac_sc_10x \
  --matrix_dir tests/data/toy_10x_mtx \
  --gtf tests/data/toy.gtf \
  --out_dir tests/tmp/readme_sc \
  --organism human \
  --genome_build hg38 \
  --groups_tsv tests/data/barcode_groups.tsv \
  --peak_summary sum_counts \
  --link_method promoter_overlap \
  --normalize l1
```
