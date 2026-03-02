# DIG Gene Set Extractors

`dig-gene-set-extractors` is the repository-level framework for building and running gene set extractors across assay types.

This repository is designed to host multiple extractor families over time.  
Current implemented family:

- `omics2geneset` (ATAC-seq, RNA-seq, and DNA methylation extractors)

## Repository Scope

This root README is intentionally framework-level:

- project structure and conventions shared across extractor families
- common output contract expectations
- where to find package-specific documentation

Detailed `omics2geneset` CLI and method behavior is split by assay:

- `docs/atac-seq2geneset.md` (ATAC practical guide)
- `docs/rna-seq-to-geneset.md` (RNA practical guide)
- `docs/scrna_cnmf_workflow.md` (scRNA cNMF preparation workflow)
- `docs/methylation2geneset.md` (DNA methylation practical guide)
- `docs/omics2geneset.md` (index/entrypoint linking assay guides)

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

ATAC-focused converters:

- `atac_bulk`
- `atac_bulk_matrix`
- `atac_sc_10x`

RNA-focused converters:

- `rna_deg`
- `rna_deg_multi`
- `rna_sc_programs`
- `sc_rna_marker`

DNA methylation converters:

- `methylation_cpg_diff`
- `methylation_dmr_regions`
- `methylation_dmr` (legacy, expects pre-aggregated gene-level inputs)

Practical guides, CLI flags, inputs, modes, and examples:

- `docs/atac-seq2geneset.md`
- `docs/rna-seq-to-geneset.md`
- `docs/scrna_cnmf_workflow.md`
- `docs/methylation2geneset.md`
- `docs/omics2geneset.md` (index)

Method notes and equations (split by assay + index):

- `docs/atac-seq_methods.tex`
- `docs/rna-seq_methods.tex`
- `docs/methylation_methods.tex`
- `docs/methods.tex` (index)

ATAC reference bundle setup:

- `docs/atac_reference_bundle.md`
- primary workflow: download one build-specific tarball (`hg19` or `hg38`), extract, and point `--resources_dir` at bundle root (no extra fetch step)
- current bundle version: `v1.1.0` split outputs: `...-atac-refdata-hg19-v1.1.0.tar.gz` and `...-atac-refdata-hg38-v1.1.0.tar.gz`
- converters auto-select build-matched resource IDs from `--genome_build`
- `omics2geneset` defaults to `--use_reference_bundle true` (opt out with `--use_reference_bundle false`)
- default ATAC preset is `connectable` and emits only high-value outputs:
  - `linked_activity` with `nearest_tss`
  - `distal_activity` with `distance_decay`
  - calibration policy `auto_prefer_ref_ubiquity_else_none` (uses reference ubiquity when available, otherwise falls back to `none`)
- `qc`, `experimental`, and `all` presets are opt-in for broader/non-recommended outputs
- when reference-backed calibration methods cannot run, converters warn and continue with available methods by default
- production guidance: run with `--resource_policy fail` after `resources status --check_schema`

Not recommended by default:

- `promoter_activity` (broad/open-chromatin bias)
- `enhancer_bias` (use via `experimental` after manual review)
- `atlas_residual` (resource-backed and score-definition sensitive; use via `experimental` or explicit `--calibration_methods`)
- `promoter_overlap` linkage outputs for primary discovery programs

## Common Output Contract

Every extractor should write:

1. `geneset.tsv`
2. `geneset.meta.json`

Optional artifacts (extractor-specific):

- `geneset.full.tsv`
- `genesets.gmt`
- `run_summary.json` / `run_summary.txt`
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
    atac-seq2geneset.md
    rna-seq-to-geneset.md
    scrna_cnmf_workflow.md
    methylation2geneset.md
    methods.tex
    atac-seq_methods.tex
    rna-seq_methods.tex
    methylation_methods.tex
    atac_reference_bundle.md
    air_gapped_install.md
  tests/
  scripts/
```

## Adding a New Extractor Family

When adding another extractor family (for example proteomics):

1. Add a new package under `src/<family_name>/`.
2. Keep extractor-specific docs under `docs/<family_name>.md`.
3. Follow the shared output contract (`geneset.tsv` + `geneset.meta.json`).
4. Add unit tests with toy fixtures under `tests/`.
5. Register a CLI entry point (either family-specific CLI or shared top-level dispatch).

## Development Notes

- Do not commit large generated artifacts under `dist/` or `data/external/`.
- For `omics2geneset` reference resources, prefer the manifest-driven resource manager workflow documented in `docs/atac-seq2geneset.md` and `docs/atac_reference_bundle.md`.
