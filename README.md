# DIG Gene Set Extractors

This repository is a framework for building and running gene set extractors across assay types.

- Repository goal: support multiple extractor families over time.
- Current implemented family: `omics2geneset`.
- First production extractors: ATAC-seq (`atac_bulk`, `atac_sc_10x`).

## Quick Start (Common Case)

```bash
python -m pip install -U pip
python -m pip install -e ".[dev]"
omics2geneset list
pytest -q
```

For offline/air-gapped setup only, see `docs/air_gapped_install.md`.

## Output Contract

All extractors write:

1. `geneset.tsv`
2. `geneset.meta.json`

ATAC extractors now default to compact program outputs:

- `geneset.tsv`: selected program (`gene_id`, `score`, `rank`, optional `weight`, optional `gene_symbol`)
- `geneset.full.tsv` (optional): full nonzero score table when `--emit_full true`
- `genesets.gmt` (optional, default on): one or more 100-500 gene sets for enrichment workflows

`geneset.meta.json` records resolved parameters, input file hashes, and summary/provenance fields.

## Program-Sized Defaults

For `atac_bulk` and `atac_sc_10x`, defaults are tuned for program extraction:

- `--select top_k`
- `--top_k 200`
- `--normalize within_set_l1`
- `--emit_full true`
- `--emit_gmt true`
- `--gmt_topk_list 200`
- `--gmt_min_genes 100`
- `--gmt_max_genes 500`

`within_set_l1` normalizes only selected genes, so `weight` sums to 1 inside the program.

## Recommended scATAC Cluster Programs

Use grouped differential extraction:

```bash
omics2geneset convert atac_sc_10x \
  --matrix_dir <10x_matrix_dir> \
  --gtf <genes.gtf> \
  --groups_tsv <barcode_to_cluster.tsv> \
  --out_dir <out_dir> \
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

Grouped runs write `group=<GROUP>/...` outputs plus `manifest.tsv`.  
`omics2geneset validate <out_dir>` validates grouped roots via the manifest.

## Choosing 100-500 Genes

- Start with `--top_k 200`.
- Use `100-150` for strict marker-style programs.
- Use `300-500` for broader pathway coverage.
- If results are too diffuse, prefer `group_vs_rest` contrast and `positive` transform before increasing `top_k`.

## About Global L1 (`--normalize l1`)

Global L1 distributes total mass across all scored genes.  
It is useful for legacy comparability, but not a probability that a gene is "active".  
For compact programs, `within_set_l1` is the recommended default.

## Documentation

- Package guide: `docs/omics2geneset.md`
- Offline install guide: `docs/air_gapped_install.md`
