import csv
import json
from pathlib import Path

import pytest

from geneset_extractors.converters import rna_deg
from geneset_extractors.core.validate import validate_output_dir


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
    padj_max = None
    pvalue_max = None
    min_abs_logfc = None
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
    gmt_topk_list = "3"
    gmt_mass_list = ""
    gmt_split_signed = True
    gmt_emit_abs = False
    gmt_source = "full"
    emit_small_gene_sets = True


def test_rna_deg_converter_end_to_end(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "rna_deg")
    result = rna_deg.run(args)
    assert result["resolved_score_mode"] == "stat"

    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
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
    provenance = json.loads((Path(args.out_dir) / "geneset.provenance.json").read_text(encoding="utf-8"))
    assert meta["converter"]["parameters"]["signature_name"] == "toy"
    assert meta["converter"]["parameters"]["score_mode"] == "stat"
    assert meta["provenance"]["path"] == "geneset.provenance.json"
    assert provenance["file_type"] == "provenance"
    assert provenance["focus_node_id"] == meta["provenance"]["focus_node_id"]
    assert any(node["kind"] == "geneset" for node in provenance["nodes"])


def test_rna_deg_provenance_overlay_injects_public_links(tmp_path: Path):
    overlay_path = tmp_path / "overlay.json"
    overlay_path.write_text(
        json.dumps(
            {
                "inputs": {
                    "role:deg_tsv": {
                        "canonical_uri": "gs://dig/example/toy_deg.tsv",
                        "download_url": "https://example.org/toy_deg.tsv",
                        "landing_page_url": "https://example.org/study",
                        "access_level": "public"
                    }
                },
                "operation": {
                    "script_url": "https://example.org/notebooks/rna_deg.py",
                    "container_image": "ghcr.io/example/dig:latest"
                }
            }
        ),
        encoding="utf-8",
    )
    args = Args()
    args.out_dir = str(tmp_path / "rna_deg_overlay")
    args.emit_gmt = False
    args.provenance_overlay_json = str(overlay_path)
    rna_deg.run(args)

    provenance = json.loads((Path(args.out_dir) / "geneset.provenance.json").read_text(encoding="utf-8"))
    deg_nodes = [node for node in provenance["nodes"] if node.get("role") == "deg_tsv"]
    assert deg_nodes
    assert deg_nodes[0]["access"]["canonical_uri"] == "gs://dig/example/toy_deg.tsv"
    assert deg_nodes[0]["access"]["download_url"] == "https://example.org/toy_deg.tsv"
    assert provenance["operations"][0]["code"]["script_url"] == "https://example.org/notebooks/rna_deg.py"
    assert provenance["operations"][0]["replay"]["container_image"] == "ghcr.io/example/dig:latest"


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


def test_rna_deg_warns_when_symbols_missing_and_ensembl_like_ids(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    no_symbol = tmp_path / "no_symbol_ensembl.tsv"
    no_symbol.write_text(
        "gene_id\tstat\nENSG00000111111\t3.0\nENSG00000122222\t-4.0\nENSG00000133333\t2.5\n",
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
    assert "gene_id appears Ensembl-like" in captured.err
    assert "gmt_require_symbol=true" in captured.err


def test_rna_deg_auto_promotes_gene_id_to_symbol_when_not_ensembl_like(tmp_path: Path):
    tsv = tmp_path / "symbols_only.tsv"
    tsv.write_text(
        "gene_id\tstat\nINS\t5.0\nABCC8\t-4.0\nPCSK1\t3.0\n",
        encoding="utf-8",
    )
    args = Args()
    args.deg_tsv = str(tsv)
    args.out_dir = str(tmp_path / "promote_symbols")
    args.gmt_require_symbol = True
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "2"
    args.emit_small_gene_sets = True
    rna_deg.run(args)

    with (Path(args.out_dir) / "geneset.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows
    assert all(str(row.get("gene_symbol", "")).strip() for row in rows)
    gmt_lines = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8").strip().splitlines()
    assert gmt_lines


def test_rna_deg_default_excludes_apply_to_gene_id_when_symbol_missing(tmp_path: Path):
    tsv = tmp_path / "gene_id_only.tsv"
    tsv.write_text(
        "gene_id\tstat\nMT-ND1\t8.0\nRPLP0\t7.0\nGENE_OK\t3.0\n",
        encoding="utf-8",
    )

    args = Args()
    args.deg_tsv = str(tsv)
    args.out_dir = str(tmp_path / "filtered_gene_id")
    args.emit_gmt = False
    rna_deg.run(args)
    with (Path(args.out_dir) / "geneset.full.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {row["gene_id"] for row in rows} == {"GENE_OK"}

    args_disable = Args()
    args_disable.deg_tsv = str(tsv)
    args_disable.out_dir = str(tmp_path / "unfiltered_gene_id")
    args_disable.emit_gmt = False
    args_disable.disable_default_excludes = True
    rna_deg.run(args_disable)
    with (Path(args_disable.out_dir) / "geneset.full.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {"MT-ND1", "RPLP0", "GENE_OK"} == {row["gene_id"] for row in rows}


def test_rna_deg_duplicate_gene_policy_sum_vs_max_abs(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    tsv = tmp_path / "dups.tsv"
    tsv.write_text(
        "gene_id\tstat\nG_DUP\t2.0\nG_DUP\t-5.0\nG_OTHER\t4.0\n",
        encoding="utf-8",
    )

    args_max = Args()
    args_max.deg_tsv = str(tsv)
    args_max.out_dir = str(tmp_path / "dup_max_abs")
    args_max.emit_gmt = False
    args_max.top_k = 1
    args_max.duplicate_gene_policy = "max_abs"
    rna_deg.run(args_max)
    captured_max = capsys.readouterr()
    assert "duplicate gene_id rows" in captured_max.err
    with (Path(args_max.out_dir) / "geneset.tsv").open("r", encoding="utf-8") as fh:
        rows_max = list(csv.DictReader(fh, delimiter="\t"))
    assert rows_max[0]["gene_id"] == "G_DUP"
    assert float(rows_max[0]["score"]) == pytest.approx(-5.0)

    args_sum = Args()
    args_sum.deg_tsv = str(tsv)
    args_sum.out_dir = str(tmp_path / "dup_sum")
    args_sum.emit_gmt = False
    args_sum.top_k = 1
    args_sum.duplicate_gene_policy = "sum"
    rna_deg.run(args_sum)
    with (Path(args_sum.out_dir) / "geneset.tsv").open("r", encoding="utf-8") as fh:
        rows_sum = list(csv.DictReader(fh, delimiter="\t"))
    assert rows_sum[0]["gene_id"] == "G_OTHER"
    assert float(rows_sum[0]["score"]) == pytest.approx(4.0)


def test_rna_deg_warns_when_biotype_allowlist_has_mostly_missing_biotypes(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
):
    args = Args()
    args.out_dir = str(tmp_path / "biotype_warn")
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "2"
    args.emit_small_gene_sets = True
    rna_deg.run(args)
    captured = capsys.readouterr()
    assert "gene_biotype values are missing" in captured.err


def test_rna_deg_default_signature_name_uses_deg_tsv_stem(tmp_path: Path):
    tsv = tmp_path / "airway_conditionA_vs_B.tsv"
    tsv.write_text("gene_id\tstat\nG1\t3.0\nG2\t-2.0\n", encoding="utf-8")
    args = Args()
    args.deg_tsv = str(tsv)
    args.signature_name = "contrast"
    args.out_dir = str(tmp_path / "default_signature")
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "2"
    args.emit_small_gene_sets = True
    rna_deg.run(args)

    meta = json.loads((Path(args.out_dir) / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta["converter"]["parameters"]["signature_name"] == "airway_conditionA_vs_B"
    gmt_text = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8")
    assert "__signature=airway_conditionA_vs_B__" in gmt_text


def test_rna_deg_gmt_emit_abs_adds_abs_ranked_sets(tmp_path: Path):
    tsv = tmp_path / "mixed.tsv"
    tsv.write_text(
        "gene_id\tstat\nA\t5.0\nB\t-7.0\nC\t2.0\n",
        encoding="utf-8",
    )
    args = Args()
    args.deg_tsv = str(tsv)
    args.out_dir = str(tmp_path / "gmt_abs")
    args.gmt_emit_abs = True
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "1"
    args.emit_small_gene_sets = True
    rna_deg.run(args)

    lines = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8").strip().splitlines()
    assert any("__abs__topk=1" in line for line in lines)
    abs_line = next(line for line in lines if "__abs__topk=1" in line)
    abs_gene = abs_line.split("\t", 1)[1]
    assert abs_gene == "B"


def test_rna_deg_warns_when_split_signed_and_no_negative_scores(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
):
    tsv = tmp_path / "all_positive.tsv"
    tsv.write_text("gene_id\tstat\nA\t5.0\nB\t3.0\nC\t1.0\n", encoding="utf-8")
    args = Args()
    args.deg_tsv = str(tsv)
    args.out_dir = str(tmp_path / "split_no_neg")
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "2"
    args.emit_small_gene_sets = True
    args.gmt_split_signed = True
    rna_deg.run(args)
    captured = capsys.readouterr()
    assert "No negative scores detected" in captured.err


def test_rna_deg_supports_signed_neglog10padj_and_avg_log2fc_autodetect(tmp_path: Path):
    tsv = tmp_path / "seurat_like.tsv"
    tsv.write_text(
        "gene_id\tavg_log2FC\tFDR\nA\t2.0\t1e-6\nB\t-1.5\t1e-4\nC\t0.2\t0.9\n",
        encoding="utf-8",
    )
    args = Args()
    args.deg_tsv = str(tsv)
    args.out_dir = str(tmp_path / "signed_padj")
    args.emit_gmt = False
    args.score_mode = "signed_neglog10padj"
    rna_deg.run(args)
    with (Path(args.out_dir) / "geneset.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows[0]["gene_id"] == "A"
    assert float(rows[0]["score"]) > 0
    assert any(float(row["score"]) < 0 for row in rows)

def test_rna_deg_supports_signed_neglog10pvalue(tmp_path: Path):
    tsv = tmp_path / "signed_pvalue.tsv"
    tsv.write_text(
        "gene_id\tlogFC\tpvalue\nA\t1.5\t1e-8\nB\t-2.0\t1e-6\nC\t0.5\t0.2\n",
        encoding="utf-8",
    )
    args = Args()
    args.deg_tsv = str(tsv)
    args.out_dir = str(tmp_path / "signed_pvalue")
    args.emit_gmt = False
    args.score_mode = "signed_neglog10pvalue"
    rna_deg.run(args)
    with (Path(args.out_dir) / "geneset.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows[0]["gene_id"] == "A"
    assert float(rows[0]["score"]) > 0
    assert any(float(row["score"]) < 0 for row in rows)


def test_rna_deg_supports_logfc_score_mode(tmp_path: Path):
    tsv = tmp_path / "logfc_mode.tsv"
    tsv.write_text(
        "gene_id\tlogFC\tpadj\nA\t1.0\t1e-6\nB\t-2.5\t1e-6\nC\t0.2\t1e-6\n",
        encoding="utf-8",
    )
    args = Args()
    args.deg_tsv = str(tsv)
    args.out_dir = str(tmp_path / "logfc_mode")
    args.emit_gmt = False
    args.score_mode = "logfc"
    rna_deg.run(args)
    with (Path(args.out_dir) / "geneset.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows[0]["gene_id"] == "B"
    assert float(rows[0]["score"]) == pytest.approx(-2.5)


def test_rna_deg_warns_for_logfc_mode_without_filters(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    tsv = tmp_path / "logfc_warn.tsv"
    tsv.write_text(
        "gene_id\tlogFC\tpadj\nA\t1.0\t1e-6\nB\t-2.5\t1e-6\nC\t0.2\t1e-6\n",
        encoding="utf-8",
    )
    args = Args()
    args.deg_tsv = str(tsv)
    args.out_dir = str(tmp_path / "logfc_warn")
    args.emit_gmt = False
    args.score_mode = "logfc"
    rna_deg.run(args)
    captured = capsys.readouterr()
    assert "score_mode=logfc is being used without" in captured.err
    summary = json.loads((Path(args.out_dir) / "run_summary.json").read_text(encoding="utf-8"))
    assert any("score_mode=logfc" in warning for warning in summary["warnings"])


def test_rna_deg_warns_when_selected_program_looks_technical(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    tsv = tmp_path / "technical.tsv"
    tsv.write_text(
        "gene_id\tstat\nMT-ND1\t10\nRPLP0\t9\nRPS3\t8\nEEF1A1\t7\nHNRNPU\t6\nGENE1\t5\nGENE2\t4\nGENE3\t3\nGENE4\t2\nGENE5\t1\n",
        encoding="utf-8",
    )
    args = Args()
    args.deg_tsv = str(tsv)
    args.out_dir = str(tmp_path / "technical_warn")
    args.emit_gmt = False
    args.disable_default_excludes = True
    args.top_k = 10
    rna_deg.run(args)
    captured = capsys.readouterr()
    assert "selected program is enriched for technical/global gene families" in captured.err
    summary = json.loads((Path(args.out_dir) / "run_summary.json").read_text(encoding="utf-8"))
    assert float(summary["technical_gene_qc"]["fraction_technical_like"]) >= 0.3
    assert any("technical/global gene families" in warning for warning in summary["warnings"])


def test_rna_deg_row_filters_apply_before_gene_aggregation(tmp_path: Path):
    tsv = tmp_path / "filtered.tsv"
    tsv.write_text(
        "gene_id\tlogFC\tpadj\tpvalue\nA\t2.0\t0.001\t0.001\nB\t1.0\t0.2\t0.2\nC\t0.1\t0.001\t0.001\n",
        encoding="utf-8",
    )
    args = Args()
    args.deg_tsv = str(tsv)
    args.out_dir = str(tmp_path / "filtered_rows")
    args.emit_gmt = False
    args.score_mode = "logfc_times_neglog10p"
    args.padj_max = 0.05
    args.min_abs_logfc = 0.5
    rna_deg.run(args)
    with (Path(args.out_dir) / "geneset.full.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert [row["gene_id"] for row in rows] == ["A"]
