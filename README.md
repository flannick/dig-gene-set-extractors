# omics2geneset

Convert omics assay outputs into a standardized gene-weight vector plus machine-readable metadata.

## Quick start

```bash
../../.venv/bin/python -m pip install -e ".[dev]"
omics2geneset list
../../.venv/bin/python -m pytest -q
```

## CLI

- `omics2geneset list`
- `omics2geneset describe <converter>`
- `omics2geneset convert <converter> [flags]`
- `omics2geneset validate <output_dir>`

Converters currently implemented:
- `atac_bulk`
- `atac_sc_10x`

Supported extraction variants:
- Peak-to-gene linking: `promoter_overlap`, `nearest_tss`, `distance_decay`
- Bulk peak transforms: `signed`, `abs`, `positive`, `negative`
- scATAC peak summaries: `sum_counts`, `mean_counts`, `frac_cells_nonzero`
- Optional grouped scATAC outputs via `--groups_tsv`
