# DIG Gene Set Extractors

This repository is a framework for building and running **gene set extractors** across assay types.  
A gene set extractor maps assay-specific inputs into a common output contract:

- `geneset.tsv` (gene weights)
- `geneset.meta.json` (provenance + parameters + summary)

## Scope

- Repository-level goal: host multiple extractor families over time.
- Current implemented family: `omics2geneset` (ATAC-seq first, extensible to additional data types).

## Current Package

- Package name: `omics2geneset`
- Source: `src/omics2geneset`
- CLI: `omics2geneset`
- Detailed package docs: `docs/omics2geneset.md`

## Quick Start (Current Package)

```bash
python -m pip install -U pip setuptools wheel
python -m pip install -e ".[dev]"
omics2geneset list
pytest -q
```

## Output Contract (All Extractors)

Every extractor should write:

1. `geneset.tsv`
2. `geneset.meta.json`

`geneset.tsv` must include at least:

- `gene_id`
- `weight`

`geneset.meta.json` must include enough metadata for reproducibility:

- schema version
- resolved converter parameters
- input file paths + hashes
- summary statistics

## Adding New Extractor Families

For a new family/package (for example a non-omics extractor module):

1. Add a package under `src/` with its own CLI integration.
2. Reuse the shared output contract.
3. Add converter specs and tests with toy data.
4. Add a package-specific guide under `docs/`.

## Documentation Layout

- Repository-level guide: `README.md` (this file)
- Current package deep-dive: `docs/omics2geneset.md`

As more extractor families are implemented, add corresponding docs in `docs/`.
