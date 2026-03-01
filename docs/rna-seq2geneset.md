# RNA-seq to Gene Sets (Practical Guide)

`omics2geneset` RNA converters consume either differential expression (DE) tables or external scRNA program loadings and emit:

- `geneset.tsv` (selected signed gene program with nonnegative weights)
- `geneset.meta.json` (provenance + parameters)
- optional `geneset.full.tsv`
- optional `genesets.gmt` (directional UP/DOWN sets)

Theory and equations: `docs/rna-seq_methods.tex`.

## Install

```bash
python -m pip install -U pip
python -m pip install -e ".[dev]"
omics2geneset list
```

## Converters

- `rna_deg`: one DE contrast per input table.
- `rna_deg_multi`: many contrasts in one long table, grouped by `--comparison_column`.
- `rna_sc_programs`: grouped scRNA gene-program extraction from external factorization loadings (generic/cNMF/scHPF outputs).
- `sc_rna_marker`: legacy single-cell count-summary converter (non-DE workflow).

## Quickstart: single contrast (`rna_deg`)

```bash
omics2geneset convert rna_deg \
  --deg_tsv path/to/deg.tsv \
  --out_dir results/rna_deg_example \
  --organism human \
  --genome_build hg38

omics2geneset validate results/rna_deg_example
```

Defaults for `rna_deg`:

- `--score_mode auto` (prefer `stat`; fallback to `logfc_times_neglog10p`)
- `--duplicate_gene_policy max_abs`
- `--signature_name contrast` resolves to `deg_tsv` filename stem
- `--select top_k --top_k 200`
- `--normalize within_set_l1`
- `--emit_gmt true --gmt_split_signed true`
- `--gmt_source full` (GMT derived from full ranked table, not only selected rows)
- `--gmt_topk_list 200 --gmt_min_genes 100 --gmt_max_genes 500`

## Quickstart: multi-contrast table (`rna_deg_multi`)

```bash
omics2geneset convert rna_deg_multi \
  --deg_tsv path/to/deg_long.tsv \
  --comparison_column comparison_id \
  --out_dir results/rna_deg_multi_example \
  --organism human \
  --genome_build hg38

omics2geneset validate results/rna_deg_multi_example
```

Grouped output layout:

- `manifest.tsv`
- `comparison=<NAME>/geneset.tsv`
- `comparison=<NAME>/geneset.meta.json`
- optional per-comparison `geneset.full.tsv` and `genesets.gmt`
- optional root `genesets.gmt` (combined)

## Quickstart: scRNA program loadings (`rna_sc_programs`)

This converter does not run NMF/cNMF/scHPF itself. It ingests precomputed gene loadings from those tools.

### Generic loadings TSV

```bash
omics2geneset convert rna_sc_programs \
  --program_loadings_tsv path/to/program_loadings.tsv \
  --loadings_format wide_genes_by_program \
  --out_dir results/rna_sc_programs \
  --organism human \
  --genome_build hg38

omics2geneset validate results/rna_sc_programs
```

### cNMF convenience input

```bash
omics2geneset convert rna_sc_programs \
  --cnmf_gene_spectra_tsv path/to/name.gene_spectra_tpm.k_20.dt_0_01.txt \
  --out_dir results/rna_sc_programs_cnmf \
  --organism human \
  --genome_build hg38
```

### scHPF convenience input

```bash
omics2geneset convert rna_sc_programs \
  --schpf_gene_scores_tsv path/to/schpf_gene_scores.tsv \
  --out_dir results/rna_sc_programs_schpf \
  --organism human \
  --genome_build hg38
```

### Output layout (`rna_sc_programs`)

- `manifest.tsv`
- `program=<ID>/geneset.tsv`
- `program=<ID>/geneset.meta.json`
- optional per-program `geneset.full.tsv` and `genesets.gmt`
- optional root `genesets.gmt` (combined across programs)

## Column mapping examples

### DESeq2-style table

```bash
omics2geneset convert rna_deg \
  --deg_tsv deseq2_results.tsv \
  --score_mode auto \
  --gene_id_column gene_id \
  --logfc_column log2FoldChange \
  --padj_column padj \
  --pvalue_column pvalue \
  --out_dir results/deseq2 \
  --organism human \
  --genome_build hg38
```

### edgeR-style table

```bash
omics2geneset convert rna_deg \
  --deg_tsv edger_results.tsv \
  --score_mode auto \
  --gene_id_column gene_id \
  --logfc_column logFC \
  --padj_column FDR \
  --pvalue_column PValue \
  --out_dir results/edger \
  --organism human \
  --genome_build hg38
```

### limma-style table

```bash
omics2geneset convert rna_deg \
  --deg_tsv limma_results.tsv \
  --score_mode auto \
  --gene_id_column gene_id \
  --logfc_column logFC \
  --pvalue_column P.Value \
  --padj_column adj.P.Val \
  --stat_column t \
  --out_dir results/limma \
  --organism human \
  --genome_build hg38
```

## Scoring modes

- `auto`: prefer `stat` if available, else `logfc_times_neglog10p`.
- `stat`: use a signed statistic column directly.
- `logfc_times_neglog10p`: `logFC * min(-log10(p), cap)` using `padj` if available, else `pvalue`.
- `custom_column`: use `--score_column` directly.

Selection and weights use `abs(score)` by default; the emitted `score` remains signed.

## scRNA program loading semantics (`rna_sc_programs`)

Per program and gene, let loading be `l_pg` from upstream factorization output.

- default `--score_transform positive`: `s_pg = max(l_pg, 0)`
- alternatives:
  - `signed`: keep signed `s_pg = l_pg`
  - `abs`: `s_pg = |l_pg|`
  - `negative`: `s_pg = max(-l_pg, 0)`

Defaults are conservative for interpretability:

- one primary positive-loading set per program
- `--select top_k --top_k 200`
- `--normalize within_set_l1`
- `--emit_gmt true --gmt_split_signed false`
- `--gmt_topk_list 200 --gmt_min_genes 100 --gmt_max_genes 500`

If loadings contain substantial negatives and `score_transform=positive`, the converter warns that negatives are being dropped and suggests signed/split mode.

Enable directional sets from signed loadings with:

```bash
--score_transform signed --gmt_split_signed true
```

Duplicate `gene_id` rows are aggregated using `--duplicate_gene_policy`:

- `max_abs` (default): keep the signed score with largest absolute magnitude
- `sum`: sum signed scores
- `mean`: mean signed score
- `last`: last row wins

## Filtering and annotation

Default symbol filters remove common technical genes:

- `^MT-`, `^RPL`, `^RPS`, `^mt-`, `^Rpl`, `^Rps`

Controls:

- add filters with repeatable `--exclude_gene_regex`
- disable defaults with `--disable_default_excludes`

Filter matching uses `gene_symbol` when present, and falls back to `gene_id` when symbols are missing.

Auto symbol handling when `gene_symbol` is missing:

- if `gene_id` is mostly non-Ensembl-like, `gene_symbol` is auto-promoted from `gene_id`
- if `gene_id` is mostly Ensembl-like, symbols are not auto-promoted; provide `--gtf` for mapping (or set `--gmt_require_symbol false`)

Optional GTF annotation:

- `--gtf <file.gtf(.gz)>`
- `--gtf_gene_id_field gene_id`
- `--gtf_source <label>`

GTF fills missing `gene_symbol`/`gene_biotype` where possible.

## GMT behavior

Default directional GMT:

- `--gmt_split_signed true` produces `__pos` and `__neg` sets
- names include sanitized `signature`, optional sanitized comparison label, and resolved `score_mode`
- default `--gmt_source full` means GMT sets are built from `geneset.full.tsv` ranking
- use `--gmt_source selected` to build GMT only from the selected `geneset.tsv` rows

Useful controls:

- `--gmt_prefer_symbol`, `--gmt_require_symbol`
- `--gmt_biotype_allowlist protein_coding`
- `--gmt_emit_abs true` to emit additional abs(score)-ranked gene sets
- `--emit_small_gene_sets true` for toy/debug runs

## Common pitfalls

- `score_mode=auto` needs either parseable `stat`, or parseable `logFC` + `p/padj`.
- If most symbols are missing and `--gmt_require_symbol true`, many genes may be dropped.
- If `--gmt_biotype_allowlist` is set but most rows have missing `gene_biotype`, the converter warns because the filter may not be informative.
- With small toy tables, lower `--gmt_min_genes` (or use `--emit_small_gene_sets true`) if GMT sets are skipped.

## Documentation map

- package index: `docs/omics2geneset.md`
- ATAC practical guide: `docs/atac-seq2geneset.md`
- RNA methods note: `docs/rna-seq_methods.tex`
- methods index: `docs/methods.tex`
