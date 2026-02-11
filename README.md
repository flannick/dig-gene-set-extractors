# omics2geneset

Convert omics assay outputs into a standardized gene-weight vector plus machine-readable metadata.

## Quick start

```bash
../.venv/bin/python -m pip install -e ".[dev]"
omics2geneset list
```

## CLI

- `omics2geneset list`
- `omics2geneset describe <converter>`
- `omics2geneset convert <converter> [flags]`
- `omics2geneset validate <output_dir>`

Converters currently implemented:
- `atac_bulk`
- `atac_sc_10x`
