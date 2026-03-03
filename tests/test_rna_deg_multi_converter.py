from pathlib import Path

import pytest

from geneset_extractors.converters import rna_deg_multi
from geneset_extractors.core.validate import validate_output_dir


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
    duplicate_gene_policy = "max_abs"
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
    gmt_emit_abs = False
    gmt_source = "full"
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

    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(out_dir, schema)


def test_rna_deg_multi_sanitizes_unsafe_comparison_label_in_gmt(tmp_path: Path):
    tsv = tmp_path / "unsafe_labels.tsv"
    tsv.write_text(
        "comparison_id\tgene_id\tstat\n"
        "\"bad label\\twith spaces\"\tG1\t5.0\n"
        "\"bad label\\twith spaces\"\tG2\t-4.0\n",
        encoding="utf-8",
    )

    args = Args()
    args.deg_tsv = str(tsv)
    args.out_dir = str(tmp_path / "unsafe_multi")
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "1"
    args.emit_small_gene_sets = True
    result = rna_deg_multi.run(args)
    assert result["n_groups"] == 1

    gmt_path = Path(args.out_dir) / "genesets.gmt"
    assert gmt_path.exists()
    lines = [line for line in gmt_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    assert lines
    for line in lines:
        assert line.count("\t") == 1
        set_name, genes = line.split("\t")
        assert set_name
        assert genes
        assert all(not ch.isspace() for ch in set_name)


def test_rna_deg_multi_default_signature_name_uses_deg_tsv_stem(tmp_path: Path):
    args = Args()
    args.signature_name = "contrast"
    args.out_dir = str(tmp_path / "multi_default_signature")
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "1"
    args.emit_small_gene_sets = True
    rna_deg_multi.run(args)

    gmt_text = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8")
    assert "__signature=toy_deg_long__" in gmt_text


def test_rna_deg_multi_biotype_missing_warning_emitted_once(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
):
    args = Args()
    args.out_dir = str(tmp_path / "multi_warn_once")
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "1"
    args.emit_small_gene_sets = True
    rna_deg_multi.run(args)
    captured = capsys.readouterr()
    assert captured.err.count("gene_biotype values are missing") == 1
