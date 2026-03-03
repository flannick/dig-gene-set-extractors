# DIG Gene Set Extractors

`dig-gene-set-extractors` is the repository-level framework for building and running gene set extractors across assay types.

This repository is designed to host multiple extractor families over time.  
Current implemented family:

- `omics2geneset` (ATAC-seq, RNA-seq, and DNA methylation extractors)
  - includes CNV segment extraction (`cnv_gene_extractor`)
  - exposed via neutral aliases `geneset-extractors`/`geneset_extractors`

## Repository Scope

This root README is intentionally framework-level:

- project structure and conventions shared across extractor families
- common output contract expectations
- where to find package-specific documentation

Detailed `omics2geneset` CLI and method behavior is split by assay:

- `docs/assays/atac/guide.md` (ATAC practical guide)
- `docs/assays/rnaseq/guide.md` (RNA practical guide)
- `docs/assays/rnaseq/scrna_cnmf_workflow.md` (scRNA cNMF preparation workflow)
- `docs/assays/methylation/guide.md` (DNA methylation practical guide)
- `docs/assays/cnv/guide.md` (CNV practical guide)
- `docs/geneset-extractors.md` (neutral index/entrypoint)
- `docs/omics2geneset.md` (compatibility index)

## Quick Start (Common Case)

```bash
python -m pip install -U pip
python -m pip install -e ".[dev]"
geneset-extractors list
pytest -q
```

For offline/air-gapped setup only, see `docs/air_gapped_install.md`.

If `pip install -e ".[dev]"` fails in DNS-restricted environments:

1. Use the offline flow in `docs/air_gapped_install.md`.
2. Use `scripts/install_offline.sh --wheelhouse <path>` for no-index installs.
3. If build isolation fails in a fresh venv, use the documented `--no-build-isolation`
   path with pre-bootstrapped local build tools (`scripts/bootstrap_build_tools.sh`).

## Implemented Extractor Family

### `omics2geneset`

CLI aliases (same implementation):

- `geneset-extractors`
- `geneset_extractors`
- `omics2geneset`

ATAC-focused converters:

- `atac_bulk`
- `atac_bulk_matrix`
- `atac_sc_10x`

RNA-focused converters:

- `rna_deg`
- `rna_deg_multi`
- `rna_sc_programs`
- `sc_rna_marker`

RNA workflow commands:

- `workflows scrna_cnmf_prepare` (downsample/split/filter + generate per-subset cNMF scripts)
- `workflows cnmf_select_k` (auto-select K from cNMF k-selection stats with reproducible heuristic)

DNA methylation converters:

- `methylation_cpg_diff`
- `methylation_dmr_regions`
- `methylation_dmr` (legacy, expects pre-aggregated gene-level inputs)

CNV converter:

- `cnv_gene_extractor`

Practical guides, CLI flags, inputs, modes, and examples:

- `docs/assays/atac/guide.md`
- `docs/assays/rnaseq/guide.md`
- `docs/assays/rnaseq/scrna_cnmf_workflow.md`
- `docs/assays/methylation/guide.md`
- `docs/assays/methylation/resources.md`
- `docs/assays/cnv/guide.md`
- `docs/geneset-extractors.md` (index)
- `docs/omics2geneset.md` (compatibility index)

Method notes and equations (split by assay + index):

- `docs/assays/atac/methods.tex`
- `docs/assays/rnaseq/methods.tex`
- `docs/assays/methylation/methods.tex`
- `docs/assays/cnv/methods.tex`
- `docs/methods.tex` (index)

ATAC reference bundle setup:

- `docs/assays/atac/reference_bundle.md`
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
    geneset_extractors/
    omics2geneset/
  docs/
    geneset-extractors.md
    omics2geneset.md
    assays/
      atac/
        guide.md
        methods.tex
        reference_bundle.md
      rnaseq/
        guide.md
        methods.tex
        scrna_cnmf_workflow.md
      methylation/
        guide.md
        methods.tex
        resources.md
      cnv/
        guide.md
        methods.tex
    methods.tex
    atac-seq2geneset.md        # compatibility alias
    rna-seq-to-geneset.md      # compatibility alias
    methylation2geneset.md     # compatibility alias
    cnv2geneset.md             # compatibility alias
    atac-seq_methods.tex       # compatibility alias
    rna-seq_methods.tex        # compatibility alias
    methylation_methods.tex    # compatibility alias
    cnv_methods.tex            # compatibility alias
    scrna_cnmf_workflow.md     # compatibility alias
    atac_reference_bundle.md   # compatibility alias
    air_gapped_install.md
  tests/
  scripts/
```

## Adding a New Extractor Family

When adding another extractor family (for example proteomics):

1. Add a new package under `src/<family_name>/`.
2. Keep extractor-specific docs under `docs/assays/<assay>/` with at least `guide.md` and `methods.tex`.
3. Follow the shared output contract (`geneset.tsv` + `geneset.meta.json`).
4. Add unit tests with toy fixtures under `tests/`.
5. Register a CLI entry point (either family-specific CLI or shared top-level dispatch).

## Development Notes

- Do not commit large generated artifacts under `dist/` or `data/external/`.
- For `omics2geneset` reference resources, prefer the manifest-driven resource manager workflow documented in `docs/assays/atac/guide.md` and `docs/assays/atac/reference_bundle.md`.
