import csv
import json
from pathlib import Path

import pytest

from omics2geneset.converters import rna_deg
from omics2geneset.core.validate import validate_output_dir


class Args:
    deg_tsv = "tests/data/toy_deg.tsv"
    out_dir = "tests/tmp/rna_deg"
    organism = "human"
    genome_build = "hg38"
    signature_name = "toy"
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
    gmt_topk_list = "3"
    gmt_mass_list = ""
    gmt_split_signed = True
    emit_small_gene_sets = True


def test_rna_deg_converter_end_to_end(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "rna_deg")
    result = rna_deg.run(args)
    assert result["resolved_score_mode"] == "stat"

    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)

    geneset = Path(args.out_dir) / "geneset.tsv"
    with geneset.open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows
    assert abs(sum(float(r["weight"]) for r in rows) - 1.0) < 1e-9

    gmt_lines = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8").strip().splitlines()
    assert any("__pos__" in line for line in gmt_lines)
    assert any("__neg__" in line for line in gmt_lines)

    meta = json.loads((Path(args.out_dir) / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta["converter"]["parameters"]["signature_name"] == "toy"
    assert meta["converter"]["parameters"]["score_mode"] == "stat"


def test_rna_deg_selection_uses_abs_score(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "rna_deg_abs")
    args.top_k = 1
    args.emit_gmt = False
    rna_deg.run(args)

    with (Path(args.out_dir) / "geneset.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["gene_id"] == "GENE_NEG"
    assert float(rows[0]["score"]) < 0


def test_rna_deg_default_symbol_filters_can_be_disabled(tmp_path: Path):
    args_default = Args()
    args_default.out_dir = str(tmp_path / "default_filter")
    args_default.top_k = 20
    args_default.emit_gmt = False
    rna_deg.run(args_default)
    with (Path(args_default.out_dir) / "geneset.full.tsv").open("r", encoding="utf-8") as fh:
        rows_default = list(csv.DictReader(fh, delimiter="\t"))
    assert "GENE_MITO" not in {r["gene_id"] for r in rows_default}
    assert "GENE_RIBO" not in {r["gene_id"] for r in rows_default}

    args_disabled = Args()
    args_disabled.out_dir = str(tmp_path / "disabled_filter")
    args_disabled.top_k = 20
    args_disabled.emit_gmt = False
    args_disabled.disable_default_excludes = True
    rna_deg.run(args_disabled)
    with (Path(args_disabled.out_dir) / "geneset.full.tsv").open("r", encoding="utf-8") as fh:
        rows_disabled = list(csv.DictReader(fh, delimiter="\t"))
    assert "GENE_MITO" in {r["gene_id"] for r in rows_disabled}
    assert "GENE_RIBO" in {r["gene_id"] for r in rows_disabled}


def test_rna_deg_auto_mode_errors_on_missing_columns(tmp_path: Path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("gene_id\tgene_symbol\nG1\tA\n", encoding="utf-8")

    args = Args()
    args.deg_tsv = str(bad)
    args.out_dir = str(tmp_path / "bad_out")
    args.emit_gmt = False
    with pytest.raises(ValueError, match="score_mode=auto"):
        rna_deg.run(args)


def test_rna_deg_warns_when_symbols_missing_and_required(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    no_symbol = tmp_path / "no_symbol.tsv"
    no_symbol.write_text(
        "gene_id\tstat\nG1\t3.0\nG2\t-4.0\nG3\t2.5\n",
        encoding="utf-8",
    )
    args = Args()
    args.deg_tsv = str(no_symbol)
    args.out_dir = str(tmp_path / "no_symbol_out")
    args.gmt_require_symbol = True
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "2"
    args.emit_small_gene_sets = True
    rna_deg.run(args)
    captured = capsys.readouterr()
    assert "no gene_symbol values found" in captured.err
