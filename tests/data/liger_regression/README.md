# LIGER regression fixture

This synthetic dataset has eight cells split across two donors, eight stable Ensembl-style gene IDs, and a matching GTF with distinct symbols. It is intentionally small and is suitable for exercising the explicit `gene_id`/`factor_id` output headers and the conversion-stage GTF annotation.

From the `dig-gene-set-extractors` repository root, prepare an executable test run:

```bash
geneset-extractors workflows scrna_liger_prepare \
  --matrix_tsv "$PWD/tests/data/liger_regression/counts.tsv" \
  --meta_tsv "$PWD/tests/data/liger_regression/metadata.tsv" \
  --meta_cell_id_column cell_id \
  --dataset_column dataset \
  --cell_type_column cell_type \
  --min_cells_per_cell_type 1 \
  --max_cells_per_bucket 8 \
  --max_cells_total 8 \
  --liger_fixed_k 2 \
  --liger_min_cells_per_dataset 1 \
  --liger_min_features 0 \
  --liger_min_umi 0 \
  --liger_max_mito 100 \
  --liger_top_n_genes 4 \
  --extractor_top_k 4 \
  --gtf "$PWD/tests/data/liger_regression/genes.gtf" \
  --gtf_source liger_regression_fixture \
  --organism human \
  --genome_build hg38 \
  --out_dir /tmp/liger_regression_run
```

Then run the generated scripts for the sole subset:

```bash
cd /tmp/liger_regression_run/subsets/regression_type
bash run_liger.sh
bash run_geneset_extractors_from_liger.sh
```

Check that `liger_out/regression_type/gene_loadings.tsv` begins with `gene_id` and contains `Factor_1` and `Factor_2`, and that converted `geneset.tsv` files include `gene_symbol` values such as `ALPHA1` and `BETA1`.

This fixture tests the TSV/metadata/GTF path. It does not test h5ad identifier preservation because this environment lacks the Python `anndata` and `h5py` packages needed to create an h5ad fixture.
