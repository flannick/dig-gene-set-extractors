from pathlib import Path

from omics2geneset.converters import rna_deg_multi
from omics2geneset.core.validate import validate_output_dir


class Args:
    deg_tsv = "tests/data/toy_deg_long.tsv"
    comparison_column = "comparison_id"
    out_dir = "tests/tmp/rna_deg_multi"
    organism = "human"
    genome_build = "hg38"
    signature_name = "batch"
    gene_id_column = "gene_id"
    gene_symbol_column = None
    stat_column = None
    logfc_column = None
    padj_column = None
    pvalue_column = None
    score_column = None
    score_mode = "auto"
    neglog10p_cap = 50.0
    neglog10p_eps = 1e-300
    exclude_gene_regex = None
    disable_default_excludes = False
    gtf = None
    gtf_gene_id_field = "gene_id"
    gtf_source = None
    select = "top_k"
    top_k = 200
    quantile = 0.01
    min_score = 0.0
    normalize = "within_set_l1"
    emit_full = True
    emit_gmt = True
    gmt_out = None
    gmt_prefer_symbol = True
    gmt_require_symbol = False
    gmt_biotype_allowlist = "protein_coding"
    gmt_min_genes = 1
    gmt_max_genes = 10
    gmt_topk_list = "2"
    gmt_mass_list = ""
    gmt_split_signed = True
    emit_small_gene_sets = True


def test_rna_deg_multi_grouped_output_and_validation(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "rna_deg_multi")
    result = rna_deg_multi.run(args)
    assert result["n_groups"] == 2

    out_dir = Path(args.out_dir)
    assert (out_dir / "manifest.tsv").exists()
    assert (out_dir / "genesets.gmt").exists()
    gmt_lines = (out_dir / "genesets.gmt").read_text(encoding="utf-8").strip().splitlines()
    assert any("__pos__" in line for line in gmt_lines)
    assert any("__neg__" in line for line in gmt_lines)

    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(out_dir, schema)
