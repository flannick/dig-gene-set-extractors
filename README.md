# omics2geneset

`omics2geneset` converts upstream omics outputs into a standardized, enrichment-ready gene set representation:

- `geneset.tsv`: one row per gene with numeric weight
- `geneset.meta.json`: machine-readable provenance and run metadata

The package is built around converter plugins so new assay-specific extractors can be added without changing the output contract.

## Quick Start

From this repo root:

```bash
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

## Output Contract

Every converter writes an output directory containing:

1. `geneset.tsv`
2. `geneset.meta.json`

`geneset.tsv` minimum columns:

- `gene_id`
- `weight`

Optional columns include `gene_symbol` and `rank`.

`geneset.meta.json` includes schema version, converter parameters, input file hashes, and summary stats (`n_input_features`, `n_genes`, `fraction_features_assigned`, etc.).

## Package Structure

- `src/omics2geneset/cli.py`: CLI entrypoint
- `src/omics2geneset/registry.py`: converter registry + spec lookup
- `src/omics2geneset/converters/`: converter implementations
- `src/omics2geneset/converters/specs/`: converter spec JSON files
- `src/omics2geneset/core/`: reusable linking/scoring/normalization/metadata code
- `src/omics2geneset/io/`: assay input parsers
- `tests/`: unit tests + toy data

## Adding New Extraction Types or Modes

### Add a new extraction type (new converter)

1. Implement a new converter module in `src/omics2geneset/converters/` with a `run(args)` function.
2. Add it to `CONVERTERS` in `src/omics2geneset/registry.py`.
3. Add CLI wiring in `src/omics2geneset/cli.py` under `convert` subcommands.
4. Add a converter spec JSON in `src/omics2geneset/converters/specs/<name>.json`.
5. Ensure output uses shared contract (`geneset.tsv` + `geneset.meta.json`).
6. Add toy input data and tests under `tests/`.

### Add a new extraction mode (existing converter)

1. Add/extend CLI flags (for example a new `--link_method` or summary method).
2. Implement mode logic inside the converter and/or shared `core/` modules.
3. Record resolved parameters in metadata (`converter.parameters`).
4. Add tests for default behavior and the new mode branch.

## Extractor: `atac_bulk`

### Purpose

Convert bulk ATAC-seq peaks plus peak-level signal into gene-level weights.

### Required Inputs

- `--peaks PATH`
  - BED-like intervals with at least `chrom`, `start`, `end`
- `--gtf PATH`
  - gene annotation GTF
- `--out_dir PATH`
- `--organism {human,mouse}`
- `--genome_build STRING`

Peak weights can come from:

- `--peaks_weight_column INT` (default `5`), or
- `--peak_weights_tsv PATH` with columns: `chrom`, `start`, `end`, `weight`

### Extraction Modes (Concept + Scientific Rationale)

Peak-to-gene linking (`--link_method`):

- `promoter_overlap`
  - assigns a peak to genes with promoter-window overlap
  - biological rationale: promoter-proximal accessibility is often strongest direct proxy for transcriptional regulation
- `nearest_tss`
  - assigns each peak to nearest transcription start site within max distance
  - rationale: nearest-gene heuristic is simple and robust when regulatory map is unknown
- `distance_decay`
  - distributes a peak across nearby genes with exponential decay by distance
  - rationale: distal cis-regulation exists; confidence generally drops with genomic distance

Peak weight transform (`--peak_weight_transform`):

- `signed`, `abs`, `positive`, `negative`

Normalization (`--normalize`):

- `l1` (default), `none`

### Outputs

In `--out_dir`:

- `geneset.tsv`
- `geneset.meta.json`

### Quickstart Command

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

Convert a 10x scATAC peak-by-cell matrix into one gene-weight vector per dataset or per user-defined group.

### Required Inputs

- `--matrix_dir PATH` containing:
  - `matrix.mtx` or `matrix.mtx.gz`
  - `barcodes.tsv` or `barcodes.tsv.gz`
  - `peaks.bed` or `features.tsv(.gz)` with peak coordinates
- `--gtf PATH`
- `--out_dir PATH`
- `--organism {human,mouse}`
- `--genome_build STRING`

Optional grouping input:

- `--groups_tsv PATH` with `barcode`, `group` columns

### Extraction Modes (Concept + Scientific Rationale)

Cell-to-group peak summarization (`--peak_summary`):

- `sum_counts`
  - total accessibility burden per peak in group
  - suitable when group size should affect signal magnitude
- `mean_counts`
  - average accessibility per cell
  - controls for group-size differences
- `frac_cells_nonzero`
  - fraction of cells with open peak
  - captures prevalence/penetrance of accessibility states

After peak summarization, the converter applies the same peak-to-gene linking framework as bulk ATAC (`promoter_overlap`, `nearest_tss`, `distance_decay`) and then normalizes gene weights.

### Outputs

Without `--groups_tsv`:

- `out_dir/geneset.tsv`
- `out_dir/geneset.meta.json`

With `--groups_tsv`:

- `out_dir/group=<GROUP>/geneset.tsv`
- `out_dir/group=<GROUP>/geneset.meta.json`
- `out_dir/manifest.tsv`

### Quickstart Command

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

## Notes

Additional converters are implemented in this repo and discoverable via `list`/`describe`. This README starts with the two core ATAC extractors above and can be extended with parallel sections for the remaining converters.
