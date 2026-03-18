# RNA-seq to Gene Sets (Practical Guide)

`geneset-extractors` RNA extractors consume either differential expression (DE) tables or external scRNA program loadings and emit:

- `geneset.tsv` (selected signed gene program with nonnegative weights)
- `geneset.meta.json` (summary metadata + provenance pointer)
- `geneset.provenance.json` (collapsed provenance graph for the emitted gene set)
- optional `geneset.full.tsv`
- optional `genesets.gmt` (directional UP/DOWN sets)

Theory and equations: `docs/assays/rnaseq/methods.tex`.

Upstream RNA preparation and inference now live under `workflows` / internal `preprocessing`, not in the extractor layer. In practice:

- `workflows rna_de_prepare`: counts + metadata -> feature preprocessing + DE -> standardized long DE table
- `workflows scrna_cnmf_prepare`: scRNA matrix + metadata -> cNMF-ready subsets/scripts
- `workflows cnmf_select_k`: cNMF k-selection helper
- `convert rna_deg`, `convert rna_deg_multi`, `convert rna_sc_programs`: already-scored assay result -> gene set

There are now two separate Harmonizome-related controls in the RNA stack:

- workflow-side `workflows rna_de_prepare --de_mode harmonizome`
  - changes how DE is fit, how samples are balanced before fitting, and how expression filtering is scoped before fitting
- extractor-side `convert rna_deg[_multi] --postprocess_mode harmonizome`
  - changes how an existing DE table is filtered/ranked into signatures

They solve different problems and can be used independently.

## Resources

RNA converters are dependency-light and do not require reference bundles.

- Required: none beyond input DE/program tables.
- Optional: `--gtf` to annotate `gene_symbol`/`gene_biotype` when those columns are absent.
- Optional workflows:
  - `workflows rna_de_prepare` generates standardized DE tables from bulk RNA-seq or scRNA-seq plus metadata.
  - `workflows scrna_cnmf_prepare` generates cNMF-ready subset matrices and shell scripts.

If you run without `--gtf`, GMT output may drop rows when `--gmt_require_symbol true` and symbols are missing.

## Install

```bash
python -m pip install -U pip
python -m pip install -e ".[dev]"
# Optional RNA helpers:
# python -m pip install -e ".[de_tools]"
# python -m pip install -e ".[scrna_tools]"
geneset-extractors list
```

## Converters

- `rna_deg`: one DE contrast per input table.
- `rna_deg_multi`: many contrasts in one long table, grouped by `--comparison_column`.
- `rna_sc_programs`: grouped scRNA gene-program extraction from external factorization loadings (generic/cNMF/scHPF outputs).
- `sc_rna_marker`: legacy single-cell count-summary converter (non-DE workflow).

Upstream RNA workflows:

- `rna_de_prepare`: bulk/scRNA pseudobulk DE staging workflow that writes `deg_long.tsv`
- `scrna_cnmf_prepare`: cNMF-oriented scRNA matrix preparation
- `cnmf_select_k`: k-selection helper

Dedicated DE workflow guide: `docs/assays/rnaseq/de_workflow.md`.

## Quickstart: single contrast (`rna_deg`)

```bash
geneset-extractors convert rna_deg \
  --deg_tsv path/to/deg.tsv \
  --out_dir results/rna_deg_example \
  --organism human \
  --genome_build hg38

geneset-extractors validate results/rna_deg_example
```

Defaults for `rna_deg`:

- `--score_mode auto` (prefer `stat`; fallback to `logfc_times_neglog10p`)
- optional explicit `--score_mode signed_neglog10padj` when you want FDR-driven ranking with direction taken from `stat` or `logFC`
- `--duplicate_gene_policy max_abs`
- `--signature_name contrast` resolves to `deg_tsv` filename stem
- `--select top_k --top_k 200`
- `--normalize within_set_l1`
- `--emit_gmt true --gmt_split_signed true`
- `--gmt_source full` (GMT derived from full ranked table, not only selected rows)
- `--gmt_topk_list 200 --gmt_min_genes 100 --gmt_max_genes 500`
- optional pre-aggregation row filters:
  - `--padj_max`
  - `--pvalue_max`
  - `--min_abs_logfc`

## Quickstart: multi-contrast table (`rna_deg_multi`)

```bash
geneset-extractors convert rna_deg_multi \
  --deg_tsv path/to/deg_long.tsv \
  --comparison_column comparison_id \
  --out_dir results/rna_deg_multi_example \
  --organism human \
  --genome_build hg38

geneset-extractors validate results/rna_deg_multi_example
```

## Post-processing modes

RNA DEG converters support two post-processing presets:

- `harmonizome` (default): tuned to reproduce the published Harmonizome RNA aging signature behavior more closely.
- `legacy`: preserves the prior RNA DEG extractor behavior.

### Main differences

`harmonizome` mode applies these RNA-specific settings:

- `--score_mode signed_neglog10padj`
- defaults `--padj_max` to `0.05` unless explicitly provided
- disables default symbol excludes such as `MT-`, `RPL`, and `RPS`
- selects genes by significance threshold rather than absolute top-k
- sets `--min_score 1.30103`, corresponding to `padj <= 0.05`
- builds GMT from `selected` rows, not the full ranked table
- emits top-250 signed GMT sets
- allows small emitted sets down to 5 genes
- clears the default protein-coding biotype allowlist

`legacy` mode preserves the previous defaults:

- `--score_mode auto` (prefer `stat`; fallback to `logfc_times_neglog10p`)
- `--select top_k --top_k 200`
- `--gmt_source full`
- `--gmt_topk_list 200 --gmt_min_genes 100 --gmt_max_genes 500`
- default technical-gene excludes remain enabled
- default protein-coding biotype allowlist remains enabled
- row filters such as `--padj_max`, `--pvalue_max`, and `--min_abs_logfc` remain user-controlled

### Ranking modes

RNA DEG converters now expose several explicit ranking intents through `--score_mode`:

- `signed_neglog10padj`
  - Direction comes from `stat` or `logFC`; rank magnitude comes from `-log10(padj)`.
  - Best default for library-style directional signatures and broad observational cohorts.
  - In the GTEx aging benchmark, this matched the published Harmonizome signatures best.

- `signed_neglog10pvalue`
  - Same signed significance ranking, but uses raw `pvalue` instead of `padj`.
  - Useful when users want pre-adjustment ordering within an explicit `padj`-filtered subset.
  - In the GTEx benchmark, it behaved almost identically to `signed_neglog10padj` at the top of the ranking.

- `stat`
  - Uses the signed model statistic directly.
  - Appropriate when the DE backend statistic is already the main object of interest.
  - In the GTEx benchmark, this tracked `signed_neglog10padj` closely when the upstream DE workflow was Harmonizome-like.

- `logfc_times_neglog10p`
  - Hybrid effect-size/significance ranking.
  - Tends to prioritize larger effects while still rewarding consistency.
  - Useful as an exploratory alternative, but it was clearly less concordant than `signed_neglog10padj` in the GTEx benchmark.

- `logfc`
  - Pure signed effect-size ranking.
  - Use when the user explicitly wants large-effect genes after applying their own row filters such as `--padj_max` or `--min_abs_logfc`.
  - Not recommended as the default for broad-tissue signature libraries: it drifted toward sparse high-effect outliers in the GTEx benchmark.

- `custom_column`
  - Uses a caller-supplied score column verbatim.
  - Use for external methods that already produced a final signed ranking.

Practical recommendation:

- Default library-style signatures: `--score_mode signed_neglog10padj`
- Raw-p ordering with the same intent: `--score_mode signed_neglog10pvalue`
- Exploratory hybrid: `--score_mode logfc_times_neglog10p`
- Effect-size-first signatures: `--score_mode logfc` with explicit row filters

### Example commands

Default Harmonizome-style run:

```bash
geneset-extractors convert rna_deg_multi \
  --deg_tsv path/to/deg_long.tsv \
  --comparison_column comparison_id \
  --out_dir results/rna_deg_multi_harmonizome \
  --organism human \
  --genome_build hg38
```

Legacy run:

```bash
geneset-extractors convert rna_deg_multi \
  --deg_tsv path/to/deg_long.tsv \
  --comparison_column comparison_id \
  --out_dir results/rna_deg_multi_legacy \
  --organism human \
  --genome_build hg38 \
  --postprocess_mode legacy
```

Ranking examples:

```bash
# Stable, significance-driven library-style ranking
geneset-extractors convert rna_deg \
  --deg_tsv path/to/deg.tsv \
  --out_dir results/rna_deg_signed_fdr \
  --score_mode signed_neglog10padj \
  --padj_max 0.05

# Effect-size-first ranking after explicit FDR filtering
geneset-extractors convert rna_deg \
  --deg_tsv path/to/deg.tsv \
  --out_dir results/rna_deg_logfc \
  --score_mode logfc \
  --padj_max 0.05 \
  --min_abs_logfc 0.25
```

Grouped output layout:

- `manifest.tsv`
- `comparison=<NAME>/geneset.tsv`
- `comparison=<NAME>/geneset.meta.json`
- `comparison=<NAME>/geneset.provenance.json`
- optional per-comparison `geneset.full.tsv` and `genesets.gmt`
- optional root `genesets.gmt` (combined)

The grouped `manifest.tsv` keeps `path` and now also includes `geneset_id`, `label`, `meta_path`, `provenance_path`, and `focus_node_id`.

## Choosing DE mode, covariates, and ranking

Use this as the practical decision table for new bulk RNA-seq datasets.

| Goal / dataset shape | Recommended DE mode | Covariates | Recommended ranking | Notes |
| --- | --- | --- | --- | --- |
| General DE inference in a controlled study | `modern` | Explicit known nuisance variables when present | `signed_neglog10padj` or `stat` | Uses all samples; best when the DE model is the primary object of interest. |
| Compact directional signature generation from heterogeneous observational cohorts | `harmonizome` | Explicit fixed effects such as `SEX`, tissue subsite, batch, RIN, ischemic/procurement timing when present | `signed_neglog10padj` | Conservative mode for broad tissues and imbalanced groups. |
| Exploratory larger-effect signatures after significance filtering | `modern` or `harmonizome`, depending on cohort | Explicit known nuisance variables | `logfc` with `--padj_max` and/or `--min_abs_logfc` | Effect-size-first mode is not a good default for library-style signatures. |
| Hybrid ranking that still rewards significance | Either | Explicit known nuisance variables | `logfc_times_neglog10p` | Useful exploratory compromise, but benchmarked worse than `signed_neglog10padj` for GTEx aging. |

How to choose covariates:

- Do not expect the tool to infer scientific covariates automatically.
- Use study-aware nuisance variables that are known or plausibly confounding for the assay and cohort.
- Common bulk RNA-seq examples: `SEX`, tissue subsite, batch, RIN/quality, ischemic or procurement timing, and other acquisition variables.

How to tell if settings are not working:

- top-ranked genes are dominated by technical/global families such as `MT-`, `RPL`, `RPS`, `EEF`, or `HNRNP`
- thousands of genes pass `padj <= 0.05` and the final top-250 set looks generic
- group sizes are strongly imbalanced in `modern` mode
- `logfc` mode is used without explicit row filters

The code now emits warnings for those objectively detectable cases, but the user still needs to validate marker/pathway coherence.

## Quickstart: prepare DE first, then extract

Use the workflow when you are starting from counts plus metadata rather than from a precomputed DE table:

```bash
geneset-extractors workflows rna_de_prepare \
  --modality bulk \
  --counts_tsv path/to/counts.tsv \
  --sample_metadata_tsv path/to/sample_metadata.tsv \
  --sample_id_column sample_id \
  --feature_id_column gene_id \
  --group_column condition \
  --comparison_mode condition_a_vs_b \
  --condition_a treated \
  --condition_b control \
  --backend auto \
  --out_dir results/rna_de_prepare \
  --organism human \
  --genome_build hg38
```

This writes `deg_long.tsv`, which can be passed directly into `rna_deg_multi`.

If you want a closer match to the published GTEx/Harmonizome aging notebook, use the explicit workflow preset instead of changing extractor settings alone:

```bash
geneset-extractors workflows rna_de_prepare \
  --modality bulk \
  --counts_tsv path/to/tissue_counts.tsv \
  --sample_metadata_tsv path/to/sample_metadata.tsv \
  --sample_id_column sample_id \
  --feature_id_column gene_id \
  --group_column age_bin \
  --comparisons_tsv path/to/gtex_age_comparisons.tsv \
  --stratify_by tissue \
  --de_mode harmonizome \
  --covariates SEX,SMTSD \
  --balance_seed 1 \
  --backend auto \
  --out_dir results/rna_de_harmonizome \
  --organism human \
  --genome_build hg38
```

This keeps the DE fit notebook-like by balancing each contrast to equal group sizes, allows explicit fixed-effect covariates for broad tissues, records the selected sample IDs in `comparison_selected_samples.tsv`, and writes pre/post balance counts plus resolved covariates to `comparison_audit.tsv`.

Use this preset when you want conservative, directional gene-set construction in heterogeneous tissues. Do not treat it as a universal RNA DE default:

- it discards samples deliberately
- it is bulk-only
- it rejects batch and repeated-measures designs
- if you omit explicit covariates, the workflow warns and records that choice in `prepare_summary.json`

If you want the workflow to call the extractor internally:

```bash
geneset-extractors workflows rna_de_prepare \
  ... \
  --run_extractor true \
  --extractor_padj_max 0.05 \
  --extractor_min_abs_logfc 0.5
```

Common workflow examples:

- bulk two-group DE: `docs/assays/rnaseq/de_workflow.md`
- scRNA donor-level pseudobulk by cell type: `docs/assays/rnaseq/de_workflow.md`
- GTEx-like age-bin contrasts within tissue: `docs/assays/rnaseq/de_workflow.md`

## Quickstart: scRNA program loadings (`rna_sc_programs`)

This converter ingests precomputed gene loadings only. It does not downsample cells, split by cell type, or run cNMF/scHPF itself.

### Recommended cNMF path in this repo

Use the workflow command to prepare cNMF-ready subsets and run scripts, then ingest the resulting gene spectra with `rna_sc_programs`:

```bash
geneset-extractors workflows scrna_cnmf_prepare \
  --matrix_tsv path/to/cell_by_gene_logcounts.tsv \
  --meta_tsv path/to/cell_metadata.tsv \
  --meta_cell_id_column cell_id \
  --cell_type_column cell_type \
  --donor_column donor_id \
  --split_by_cell_type true \
  --max_cells_per_bucket 200 \
  --max_cells_total 20000 \
  --cnmf_k_list auto \
  --cnmf_k auto \
  --out_dir results/scrna_cnmf_prepare

# Then run each generated subset script:
#   results/scrna_cnmf_prepare/subsets/<subset>/run_cnmf.sh
#   results/scrna_cnmf_prepare/subsets/<subset>/run_cnmf_consensus_auto_k.sh
#   results/scrna_cnmf_prepare/subsets/<subset>/run_geneset_extractors_from_cnmf.sh

# After cNMF consensus, ingest gene spectra:
geneset-extractors convert rna_sc_programs \
  --cnmf_gene_spectra_tsv path/to/<name>.gene_spectra_tpm.k_<K>.dt_<...>.txt \
  --out_dir results/rna_sc_programs_cnmf \
  --organism human \
  --genome_build hg38
```

Detailed workflow page: `docs/assays/rnaseq/scrna_cnmf_workflow.md`.

If your matrix is gene x cell, use:

```bash
--matrix_orientation gene_by_cell --matrix_gene_id_column gene_id
```

Override examples:

```bash
# Fixed K grid and fixed consensus K
geneset-extractors workflows scrna_cnmf_prepare \
  ... \
  --cnmf_k_list "10 15 20 25" \
  --cnmf_k 20
```

### Generic loadings TSV

```bash
geneset-extractors convert rna_sc_programs \
  --program_loadings_tsv path/to/program_loadings.tsv \
  --loadings_format wide_genes_by_program \
  --out_dir results/rna_sc_programs \
  --organism human \
  --genome_build hg38

geneset-extractors validate results/rna_sc_programs
```

### cNMF convenience input

```bash
geneset-extractors convert rna_sc_programs \
  --cnmf_gene_spectra_tsv path/to/name.gene_spectra_tpm.k_20.dt_0_01.txt \
  --out_dir results/rna_sc_programs_cnmf \
  --organism human \
  --genome_build hg38
```

### scHPF convenience input

```bash
geneset-extractors convert rna_sc_programs \
  --schpf_gene_scores_tsv path/to/schpf_gene_scores.tsv \
  --out_dir results/rna_sc_programs_schpf \
  --organism human \
  --genome_build hg38
```

### Output layout (`rna_sc_programs`)

- `manifest.tsv`
- `program=<ID>/geneset.tsv`
- `program=<ID>/geneset.meta.json`
- `program=<ID>/geneset.provenance.json`
- optional per-program `geneset.full.tsv` and `genesets.gmt`
- optional root `genesets.gmt` (combined across programs)

## Optional provenance overlay

If your DE table or program loadings live behind public object storage or a portal landing page, attach that information without changing the extractor inputs:

```bash
geneset-extractors convert rna_deg \
  --deg_tsv path/to/deg.tsv \
  --out_dir results/rna_deg_example \
  --organism human \
  --genome_build hg38 \
  --provenance_overlay_json provenance_overlay.json
```

Typical overlay entries map an input path or `role:<input_role>` to `canonical_uri`, `download_url`, `landing_page_url`, or `persistent_id`, and can also set operation-level `script_url`, `notebook_url`, `container_image`, or `workspace_template_url`.

Local-only provenance is still valid. Public URLs are only emitted when the converter or overlay explicitly knows them.

## Best practices for large scRNA maps (recommended upstream factorization workflow)

When working with large, heterogeneous single-cell atlases, the quality of inferred programs usually depends more on upstream factorization strategy than on downstream conversion.

Recommended upstream workflow:

- Fit programs within cell types or a small set of broad compartments rather than one global run across all cell states.
- Downsample before factorization:
  - cap total cells
  - cap cells per donor (and per cell type if available)
- If donor metadata is present, avoid donor imbalance so one donor does not dominate factorization.
- Keep gene filtering policy stable between factorization and conversion (especially mitochondrial/ribosomal handling).
- Start with moderate `K` values and increase only if programs remain interpretable.
- cNMF expects nonnegative expression matrices and should not include zero-total cells or zero-total genes. The `workflows scrna_cnmf_prepare` command enforces these filters.

Large atlas best-practice checklist:

- Split by cell type when the atlas is heterogeneous.
- Cap cells per `(cell_type, donor)` bucket to avoid donor dominance.
- Keep `--seed` fixed for reproducibility.
- Prefer counts for cNMF when available; if using logcounts, keep values nonnegative and use the generated `--densify` script option as needed.

Supported loadings artifacts for `rna_sc_programs`:

- cNMF gene spectra (`gene_spectra_tpm` or `gene_spectra_score`)
- scHPF gene scores table
- Generic TSV in either wide genes-by-program or long tidy format

Minimal long-tidy TSV example:

```tsv
program_id	gene_id	loading
P1	INS	2.41
P1	MAFA	1.73
```

If you see warnings about very large numbers of programs or parsed values, split conversion by cell type and/or reduce factor count before rerunning.

## Column mapping examples

### DESeq2-style table

```bash
geneset-extractors convert rna_deg \
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
geneset-extractors convert rna_deg \
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
geneset-extractors convert rna_deg \
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

Repeated `(program_id, gene_id)` rows are aggregated additively at parse time.

## Optional helper script for quick NMF loadings

For exploratory workflows (outside cNMF), the repo includes an optional convenience script:

```bash
python scripts/scrna_quick_nmf.py --help
```

It reads:

- a cell x gene expression TSV
- optional metadata TSV (`cell_id` plus optional `cell_type` and `donor_id`)

Then it:

- downsamples per `(cell_type, donor_id)` bucket
- fits sklearn NMF
- writes a wide genes-by-program table compatible with `rna_sc_programs`

Optional dependencies (not required for base `geneset-extractors`):

```bash
python -m pip install -e ".[scrna_tools]"
```

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

- package index: `docs/geneset-extractors.md`
- ATAC practical guide: `docs/assays/atac/guide.md`
- scRNA cNMF workflow: `docs/assays/rnaseq/scrna_cnmf_workflow.md`
- RNA methods note: `docs/assays/rnaseq/methods.tex`
- methods index: `docs/methods.tex`
