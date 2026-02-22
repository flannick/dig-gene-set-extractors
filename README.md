# DIG Gene Set Extractors

This repository is a framework for building and running gene set extractors across assay types.

- Repository goal: support multiple extractor families over time.
- Current implemented family: `omics2geneset`.
- First production extractors: ATAC-seq (`atac_bulk`, `atac_bulk_matrix`, `atac_sc_10x`).

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
- `--link_method all`
- `--program_preset default`
- `--gmt_topk_list 100,200,500`
- `--gmt_mass_list 0.5,0.8,0.9`
- `--gmt_require_symbol true`
- `--gmt_biotype_allowlist protein_coding`
- `--gmt_min_genes 100`
- `--gmt_max_genes 500`

`within_set_l1` normalizes only selected genes, so `weight` sums to 1 inside the program.
`program_preset=default` emits multiple ATAC program families in GMT:

- bulk: `linked_activity`, `promoter_activity`, `distal_activity`, `enhancer_bias`
- scATAC: same plus `tfidf_distal`

Optional resource-backed methods:

- `ref_ubiquity_penalty`
- `atlas_residual`

Enable them with `--program_preset all` or explicit `--program_methods ...`.
Default preset runs (`--program_preset default`) do not require any external resource downloads.

`program_preset=connectable` prioritizes direction-aware programs for semantic handles (condition/phenotype/time/genotype):

- bulk/sc: keeps `linked_activity`, `distal_activity`, `enhancer_bias` (promoter-only program is excluded by default).
- scATAC with `--cell_metadata_tsv` + `--condition_column`: auto-switches to `condition_within_group` and emits OPEN/CLOSE programs.
- scATAC without condition metadata: falls back to `group_vs_rest`.

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
  --program_preset default \
  --peak_weight_transform positive \
  --select top_k \
  --top_k 200 \
  --normalize within_set_l1
```

## Connectable Contrast Programs

### scATAC: condition-within-group OPEN/CLOSE

Use this when you have cell metadata with condition labels:

```bash
omics2geneset convert atac_sc_10x \
  --matrix_dir <10x_matrix_dir> \
  --gtf <genes.gtf> \
  --groups_tsv <barcode_to_cluster.tsv> \
  --cell_metadata_tsv <cell_metadata.tsv> \
  --condition_column condition \
  --condition_a treated \
  --condition_b control \
  --min_cells_per_condition 50 \
  --out_dir <out_dir> \
  --organism human \
  --genome_build hg38 \
  --program_preset connectable
```

Expected outputs:

- `group=<GROUP>/geneset.tsv` (selected direction; OPEN by default, CLOSE if `--peak_weight_transform negative`)
- `group=<GROUP>/genesets.gmt` with both `direction=OPEN` and `direction=CLOSE` sets
- grouped root `genesets.gmt` and `manifest.tsv`

### bulk ATAC across samples: OPEN/CLOSE differential programs

Use `atac_bulk_matrix` when you have a peak-by-sample matrix plus sample metadata:

```bash
omics2geneset convert atac_bulk_matrix \
  --peak_matrix_tsv <peak_by_sample.tsv> \
  --peak_bed <peaks.bed> \
  --sample_metadata_tsv <sample_metadata.tsv> \
  --sample_id_column sample_id \
  --condition_column condition \
  --condition_a case \
  --condition_b control \
  --gtf <genes.gtf> \
  --out_dir <out_dir> \
  --organism human \
  --genome_build hg38 \
  --program_preset connectable
```

Expected outputs:

- `geneset.tsv` (selected direction)
- `genesets.gmt` with both OPEN and CLOSE sets and contrast labels in set names
- `geneset.meta.json` with `condition_a/b`, `n_samples_a/b`, contrast metric and pseudocount

Grouped runs write `group=<GROUP>/...` outputs plus `manifest.tsv`.  
`omics2geneset validate <out_dir>` validates grouped roots via the manifest.

External linkage (ABC/Hi-C style) is supported via:

- `--link_method external` or combined modes such as `--link_method all,external`
- `--region_gene_links_tsv <path>`

Standardized external linkage TSV format:

- required header columns: `chrom`, `start`, `end`, `gene_id`, `link_weight`
- coordinates are 0-based half-open intervals
- `gene_id` should match GTF gene IDs (versioned or unversioned Ensembl IDs are both accepted)

Optional resource-backed program flags:

- `--resource_policy {skip,fail}` (default `skip`)
- `--resources_manifest <manifest.json>` (optional override)
- `--resources_dir <resource_cache_dir>` (optional override)
- `--ref_ubiquity_resource_id <id>` and `--atlas_resource_id <id>`
- `--atlas_metric {logratio,zscore}` and `--atlas_eps <float>`

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
- ATAC reference bundle guide: `docs/atac_reference_bundle.md`

## Optional Resource Catalog

For methods that need external references (for example cCRE ubiquity references), use:

```bash
omics2geneset resources list
omics2geneset resources status --fast
omics2geneset resources describe ccre_ubiquity_hg38
omics2geneset resources fetch --preset atac_default_optional
omics2geneset resources manifest-validate
```

Set `OMICS2GENESET_RESOURCES_DIR` to override the default cache location.
By default, `resources fetch` skips entries without a URL and reports them as `manual`.

Use `--manifest <path>` to provide a custom manifest.
Default mode is overlay merge with the bundled manifest; set `--manifest_mode replace` to use only your file.

Resource file formats used by resource-backed program methods:

- `ref_ubiquity_penalty` table:
  - required columns: `chrom`, `start`, `end`, plus either:
    - `idf_ref`, or
    - `df_ref` and `n_ref` (idf computed as `log((n_ref+1)/(df_ref+1))+1`)
- `atlas_residual` table:
  - required columns: `gene_id`, `median_score`, `mad_score`
  - accepted aliases: `median` / `mad`

Resource IDs currently referenced by ATAC methods:

- `ccre_ubiquity_hg38`, `ccre_ubiquity_mm10`: used by `ref_ubiquity_penalty`
- `atac_reference_profiles_hg38`, `atac_reference_profiles_mm10`: used by `atlas_residual`

Each resource entry tracks `provider`, `stable_id`, `version`, `license`, `filename`, optional `url`, and optional `sha256`.

## ATAC Reference Bundle (Download + Local Usage)

If you have the published bundle URL, see `docs/atac_reference_bundle.md` for the full workflow.

If you built the bundle locally and do not have web hosting yet, you can use Unix paths directly:

```bash
# 1) Create a local manifest with absolute Unix paths in each resource "url" field.
python scripts/make_local_resources_manifest.py \
  --bundle-root /path/to/bundle \
  --out /tmp/omics2geneset.local_resources.json

# 2) Copy resources from those local paths into your resource cache.
omics2geneset resources fetch \
  --manifest /tmp/omics2geneset.local_resources.json \
  --manifest_mode replace \
  --preset atac_default_optional
```

Notes:

- `resources fetch` accepts local Unix filesystem paths in `url` (no `file://` URI required).
- You can then run converters with `--resources_manifest /tmp/omics2geneset.local_resources.json`.
