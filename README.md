# DIG Gene Set Extractors

`dig-gene-set-extractors` is the repository-level framework for building and running gene set extractors across assay types.

This repository is designed to host multiple extractor families over time.  
Current implemented family:

- `omics2geneset` (ATAC-seq extractors)

## Repository Scope

This root README is intentionally framework-level:

- project structure and conventions shared across extractor families
- common output contract expectations
- where to find package-specific documentation

Detailed `omics2geneset` CLI and method behavior lives in `docs/omics2geneset.md`.

## Quick Start (Common Case)

```bash
python -m pip install -U pip
python -m pip install -e ".[dev]"
omics2geneset list
pytest -q
```

For offline/air-gapped setup only, see `docs/air_gapped_install.md`.

## Implemented Extractor Family

### `omics2geneset`

ATAC-focused converters currently available:

- `atac_bulk`
- `atac_bulk_matrix`
- `atac_sc_10x`

Full package guide, CLI flags, inputs, modes, and examples:

- `docs/omics2geneset.md`

Method notes and equations:

- `docs/methods.tex`

ATAC reference bundle setup:

- `docs/atac_reference_bundle.md`

## Common Output Contract

Every extractor should write:

1. `geneset.tsv`
2. `geneset.meta.json`

Optional artifacts (extractor-specific):

- `geneset.full.tsv`
- `genesets.gmt`
- grouped outputs with `manifest.tsv`

At minimum, output metadata should record:

- extractor/converter identity and resolved algorithmic parameters
- input file provenance (including hashes)
- summary stats needed for downstream reproducibility checks

## Repository Layout

```text
dig-gene-set-extractors/
  README.md
  pyproject.toml
  src/
    omics2geneset/
  docs/
    omics2geneset.md
    methods.tex
    atac_reference_bundle.md
    air_gapped_install.md
  tests/
  scripts/
```

## Adding a New Extractor Family

When adding another extractor family (for example RNA-seq or proteomics):

1. Add a new package under `src/<family_name>/`.
2. Keep extractor-specific docs under `docs/<family_name>.md`.
3. Follow the shared output contract (`geneset.tsv` + `geneset.meta.json`).
4. Add unit tests with toy fixtures under `tests/`.
5. Register a CLI entry point (either family-specific CLI or shared top-level dispatch).

## Development Notes

- Do not commit large generated artifacts under `dist/` or `data/external/`.
- For `omics2geneset` reference resources, prefer the manifest-driven resource manager workflow documented in `docs/omics2geneset.md` and `docs/atac_reference_bundle.md`.
