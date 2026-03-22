import csv
import json
from pathlib import Path

import pytest

from geneset_extractors.converters import drug_response_screen
from geneset_extractors.core.validate import validate_output_dir
from tests.provenance_helpers import (
    assert_node_has_structured_resource_metadata,
    file_node_for_role,
    load_provenance,
)


class Args:
    out_dir = "tests/tmp/drug_response"
    organism = "human"
    genome_build = "hg38"
    dataset_label = "toy_drug_response"
    program_preset = "connectable"

    response_tsv = "tests/data/toy_drug_response.tsv"
    response_delimiter = "\t"
    sample_id_column = "sample_id"
    drug_id_column = "drug_id"
    response_column = "response"

    drug_targets_tsv = "tests/data/toy_drug_targets.tsv"
    targets_delimiter = "\t"
    targets_drug_id_column = "drug_id"
    targets_gene_symbol_column = "gene_symbol"
    targets_weight_column = "weight"
    targets_source_column = "source"

    sample_metadata_tsv = None
    sample_metadata_delimiter = "\t"
    sample_metadata_sample_id_column = "sample_id"
    group_column = "group"

    groups_tsv = "tests/data/toy_drug_groups.tsv"
    groups_delimiter = "\t"
    groups_sample_id_column = "sample_id"
    groups_group_column = "group"

    case_control_tsv = None
    case_control_delimiter = "\t"
    case_control_sample_id_column = "sample_id"
    case_control_is_case_column = "is_case"
    case_control_within_group = False

    prism_matrix_csv = None
    prism_treatment_info_csv = None
    prism_cell_line_info_csv = None
    prism_matrix_sample_id_column = "row_name"
    prism_cell_line_sample_id_column = "row_name"
    prism_column_name_column = "column_name"
    prism_broad_id_column = "broad_id"
    prism_target_column = "target"

    response_metric = "logfold_change"
    response_direction = None
    response_transform = "robust_z_mad"
    contrast_method = "group_vs_rest"
    min_group_size = 1
    scoring_model = "target_weighted_sum"
    sparse_alpha = 0.01
    response_ubiquity_penalty = None
    ubiquity_penalty = "fraction_active"
    target_ubiquity_penalty = "idf"
    ubiquity_tau = 1.0
    ubiquity_epsilon = 0.05
    polypharm_downweight = True
    polypharm_t0 = 5
    max_targets_per_drug = 50
    target_promiscuity_policy = "warn"
    max_programs = 50

    select = "top_k"
    top_k = 2
    quantile = 0.01
    min_score = 0.0
    normalize = "within_set_l1"
    emit_full = True

    target_aliases_tsv = None
    alias_delimiter = "\t"
    alias_column = "alias"
    alias_gene_symbol_column = "gene_symbol"
    drug_blacklist_tsv = None
    drug_blacklist_delimiter = "\t"
    resources_manifest = None
    resources_dir = None
    resource_policy = "skip"
    target_aliases_resource_id = None
    drug_blacklist_resource_id = None
    drug_alias_map_resource_id = None
    target_edges_resource_id = None
    target_ubiquity_resource_id = None
    compound_qc_resource_id = None
    use_compound_qc_bundle = None
    include_blacklisted_compounds = False

    gtf = None
    gtf_source = None
    gtf_gene_id_field = "gene_id"

    emit_gmt = True
    gmt_out = None
    gmt_prefer_symbol = True
    gmt_require_symbol = False
    gmt_biotype_allowlist = "protein_coding"
    gmt_min_genes = 1
    gmt_max_genes = 10
    gmt_topk_list = "2"
    gmt_mass_list = ""
    gmt_split_signed = None
    gmt_format = "dig2col"
    emit_small_gene_sets = True
    provenance_overlay_json = None


def _read_manifest_rows(out_dir: Path) -> list[dict[str, str]]:
    with (out_dir / "manifest.tsv").open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _read_geneset_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_drug_response_generic_group_vs_rest_signed_gmt(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "drug_response_generic")
    result = drug_response_screen.run(args)
    assert result["n_groups"] == 2

    out_dir = Path(args.out_dir)
    root_summary = json.loads((out_dir / "run_summary.json").read_text(encoding="utf-8"))
    assert int(root_summary["n_response_rows"]) > 0
    assert int(root_summary["n_cell_lines"]) > 0
    assert int(root_summary["n_compounds"]) > 0
    root_summary_txt = (out_dir / "run_summary.txt").read_text(encoding="utf-8")
    assert "n_response_rows: 0" not in root_summary_txt
    rows = _read_manifest_rows(out_dir)
    assert len(rows) == 2
    for row in rows:
        p = out_dir / str(row["path"])
        assert (p / "geneset.tsv").exists()
        gmt_path = p / "genesets.gmt"
        assert gmt_path.exists()
        text = gmt_path.read_text(encoding="utf-8")
        assert "__pos__" in text
        assert "__neg__" in text

    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(out_dir, schema)


def test_drug_response_bundle_resources_are_enriched_in_provenance(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "drug_response_bundle")
    args.drug_targets_tsv = None
    args.resources_dir = "tests/data/drug_response_bundle"
    result = drug_response_screen.run(args)
    assert result["n_groups"] >= 1

    manifest_rows = _read_manifest_rows(Path(args.out_dir))
    provenance = load_provenance(Path(args.out_dir) / manifest_rows[0]["path"])
    node = file_node_for_role(provenance, "target_ubiquity_bundle_tsv")
    assert_node_has_structured_resource_metadata(node)


def test_drug_response_prism_convenience_mode(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "drug_response_prism")
    args.response_tsv = None
    args.drug_targets_tsv = None
    args.groups_tsv = None
    args.contrast_method = None
    args.prism_matrix_csv = "tests/data/toy_prism_matrix.csv"
    args.prism_treatment_info_csv = "tests/data/toy_prism_treatment_info.csv"
    args.prism_cell_line_info_csv = "tests/data/toy_prism_cell_line_info.csv"
    result = drug_response_screen.run(args)
    assert result["n_groups"] == 2
    out_dir = Path(args.out_dir)
    assert (out_dir / "manifest.tsv").exists()
    assert (out_dir / "genesets.gmt").exists()


def test_drug_response_ubiquity_penalty_changes_top_gene(tmp_path: Path):
    args_none = Args()
    args_none.out_dir = str(tmp_path / "drug_response_no_penalty")
    args_none.response_tsv = "tests/data/toy_drug_response_ubiquity.tsv"
    args_none.drug_targets_tsv = "tests/data/toy_drug_targets_ubiquity.tsv"
    args_none.groups_tsv = None
    args_none.contrast_method = "none"
    args_none.response_metric = "other"
    args_none.response_direction = "higher_is_more_sensitive"
    args_none.response_transform = "none"
    args_none.ubiquity_penalty = "none"
    args_none.max_programs = 1
    args_none.top_k = 1
    args_none.gmt_topk_list = "1"
    drug_response_screen.run(args_none)

    args_pen = Args()
    args_pen.out_dir = str(tmp_path / "drug_response_with_penalty")
    args_pen.response_tsv = "tests/data/toy_drug_response_ubiquity.tsv"
    args_pen.drug_targets_tsv = "tests/data/toy_drug_targets_ubiquity.tsv"
    args_pen.groups_tsv = None
    args_pen.contrast_method = "none"
    args_pen.response_metric = "other"
    args_pen.response_direction = "higher_is_more_sensitive"
    args_pen.response_transform = "none"
    args_pen.ubiquity_penalty = "fraction_active"
    args_pen.max_programs = 1
    args_pen.top_k = 1
    args_pen.gmt_topk_list = "1"
    drug_response_screen.run(args_pen)

    out_none = Path(args_none.out_dir)
    out_pen = Path(args_pen.out_dir)
    row_none = _read_manifest_rows(out_none)[0]
    row_pen = _read_manifest_rows(out_pen)[0]
    gene_none = _read_geneset_rows(out_none / row_none["path"] / "geneset.tsv")[0]["gene_id"]
    gene_pen = _read_geneset_rows(out_pen / row_pen["path"] / "geneset.tsv")[0]["gene_id"]
    assert gene_none == "GENEUB"
    assert gene_pen == "GENESEL"


def test_drug_response_low_target_coverage_warning(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    args = Args()
    args.out_dir = str(tmp_path / "drug_response_low_coverage")
    args.response_tsv = "tests/data/toy_drug_response.tsv"
    args.drug_targets_tsv = "tests/data/toy_drug_targets_low_coverage.tsv"
    drug_response_screen.run(args)
    captured = capsys.readouterr()
    assert "low drug-to-target coverage" in captured.err


def _write_tiny_group_fixture(base: Path) -> tuple[Path, Path, Path]:
    response = base / "response.tsv"
    groups = base / "groups.tsv"
    targets = base / "targets.tsv"
    response.write_text(
        "\n".join(
            [
                "sample_id\tdrug_id\tresponse",
                "S1\tD1\t-2.0",
                "S1\tD2\t-1.0",
                "S1\tD3\t-0.5",
                "S2\tD1\t-1.5",
                "S2\tD2\t-0.8",
                "S2\tD3\t-0.4",
                "S3\tD1\t-1.4",
                "S3\tD2\t-0.7",
                "S3\tD3\t-0.3",
                "S4\tD1\t-1.2",
                "S4\tD2\t-0.6",
                "S4\tD3\t-0.2",
                "S5\tD1\t-0.2",
                "S5\tD2\t-2.1",
                "S5\tD3\t-0.7",
                "S6\tD1\t-0.1",
                "S6\tD2\t-2.0",
                "S6\tD3\t-0.6",
                "S7\tD1\t-0.2",
                "S7\tD2\t-1.9",
                "S7\tD3\t-0.5",
                "S8\tD1\t-0.3",
                "S8\tD2\t-1.8",
                "S8\tD3\t-0.4",
                "S9\tD1\t-0.2",
                "S9\tD2\t-1.7",
                "S9\tD3\t-0.3",
                "S10\tD1\t-0.1",
                "S10\tD2\t-1.6",
                "S10\tD3\t-0.2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    groups.write_text(
        "\n".join(
            [
                "sample_id\tgroup",
                "S1\tG1",
                "S2\tG2",
                "S3\tG2",
                "S4\tG2",
                "S5\tG3",
                "S6\tG3",
                "S7\tG3",
                "S8\tG3",
                "S9\tG3",
                "S10\tG3",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    target_lines = ["drug_id\tgene_symbol\tweight", "D1\tGENEA\t1.0", "D2\tGENEB\t1.0", "D3\tGENEC\t1.0", "D4\t\t1.0"]
    for idx in range(1, 101):
        target_lines.append(f"D5\tPROM{idx:03d}\t1.0")
    targets.write_text("\n".join(target_lines) + "\n", encoding="utf-8")
    return response, groups, targets


def test_drug_response_default_min_group_size_skips_tiny_groups(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    response, groups, targets = _write_tiny_group_fixture(tmp_path)
    args = Args()
    args.out_dir = str(tmp_path / "drug_response_tiny_groups")
    args.response_tsv = str(response)
    args.groups_tsv = str(groups)
    args.drug_targets_tsv = str(targets)
    args.min_group_size = 5
    args.target_promiscuity_policy = "warn"
    args.max_targets_per_drug = 50
    args.top_k = 3
    args.gmt_topk_list = "3"
    args.gmt_min_genes = 1
    args.gmt_max_genes = 20
    result = drug_response_screen.run(args)
    assert result["n_groups"] == 1
    out_dir = Path(args.out_dir)
    manifest_rows = _read_manifest_rows(out_dir)
    assert len(manifest_rows) == 1
    assert manifest_rows[0]["group"] == "G3"
    group_qc_rows = list(csv.DictReader((out_dir / "group_qc.tsv").open("r", encoding="utf-8"), delimiter="\t"))
    skipped = [r for r in group_qc_rows if r["reason_if_skipped"] == "min_group_size"]
    assert len(skipped) == 2
    assert {r["group"] for r in skipped} == {"G1", "G2"}
    run_summary = json.loads((out_dir / "run_summary.json").read_text(encoding="utf-8"))
    assert int(run_summary["n_groups_skipped"]) == 2
    captured = capsys.readouterr()
    assert "group_skipped group=G1 n=1 reason=min_group_size" in captured.err
    assert "promiscuous target list for drug=D5" in captured.err


def test_drug_response_gmt_classic_format(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "drug_response_classic_gmt")
    args.gmt_format = "classic"
    result = drug_response_screen.run(args)
    assert result["n_groups"] == 2
    out_dir = Path(args.out_dir)
    row = _read_manifest_rows(out_dir)[0]
    gmt_path = out_dir / row["path"] / "genesets.gmt"
    line = gmt_path.read_text(encoding="utf-8").strip().splitlines()[0]
    assert len(line.split("\t")) >= 3


def test_target_ubiquity_penalty_idf_changes_weights(tmp_path: Path):
    response = tmp_path / "resp.tsv"
    targets = tmp_path / "targets.tsv"
    response.write_text(
        "\n".join(
            [
                "sample_id\tdrug_id\tresponse",
                "S1\tD1\t0.7",
                "S1\tD2\t0.7",
                "S1\tD3\t0.6",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    targets.write_text(
        "\n".join(
            [
                "drug_id\tgene_symbol\tweight",
                "D1\tGENEUB\t1.0",
                "D2\tGENEUB\t1.0",
                "D3\tGENESEL\t1.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    args_none = Args()
    args_none.out_dir = str(tmp_path / "target_ubiq_none")
    args_none.response_tsv = str(response)
    args_none.drug_targets_tsv = str(targets)
    args_none.groups_tsv = None
    args_none.contrast_method = "none"
    args_none.response_metric = "other"
    args_none.response_direction = "higher_is_more_sensitive"
    args_none.response_transform = "none"
    args_none.response_ubiquity_penalty = "none"
    args_none.ubiquity_penalty = "none"
    args_none.target_ubiquity_penalty = "none"
    args_none.max_programs = 1
    args_none.top_k = 1
    args_none.gmt_topk_list = "1"
    drug_response_screen.run(args_none)

    args_idf = Args()
    args_idf.out_dir = str(tmp_path / "target_ubiq_idf")
    args_idf.response_tsv = str(response)
    args_idf.drug_targets_tsv = str(targets)
    args_idf.groups_tsv = None
    args_idf.contrast_method = "none"
    args_idf.response_metric = "other"
    args_idf.response_direction = "higher_is_more_sensitive"
    args_idf.response_transform = "none"
    args_idf.response_ubiquity_penalty = "none"
    args_idf.ubiquity_penalty = "none"
    args_idf.target_ubiquity_penalty = "idf"
    args_idf.max_programs = 1
    args_idf.top_k = 1
    args_idf.gmt_topk_list = "1"
    drug_response_screen.run(args_idf)

    row_none = _read_manifest_rows(Path(args_none.out_dir))[0]
    row_idf = _read_manifest_rows(Path(args_idf.out_dir))[0]
    gene_none = _read_geneset_rows(Path(args_none.out_dir) / row_none["path"] / "geneset.tsv")[0]["gene_id"]
    gene_idf = _read_geneset_rows(Path(args_idf.out_dir) / row_idf["path"] / "geneset.tsv")[0]["gene_id"]
    assert gene_none == "GENEUB"
    assert gene_idf == "GENESEL"


def test_broad_pharmacology_warning_emitted(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    response = tmp_path / "resp.tsv"
    targets = tmp_path / "targets.tsv"
    response.write_text(
        "\n".join(
            [
                "sample_id\tdrug_id\tresponse",
                "S1\tD1\t1.0",
                "S1\tD2\t0.9",
                "S1\tD3\t0.8",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    targets.write_text(
        "\n".join(
            [
                "drug_id\tgene_symbol\tweight",
                "D1\tGPR12\t1.0",
                "D2\tHTR2A\t1.0",
                "D3\tDRD2\t1.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    args = Args()
    args.out_dir = str(tmp_path / "broad_warning")
    args.response_tsv = str(response)
    args.drug_targets_tsv = str(targets)
    args.groups_tsv = None
    args.contrast_method = "none"
    args.response_metric = "other"
    args.response_direction = "higher_is_more_sensitive"
    args.response_transform = "none"
    args.max_programs = 1
    args.top_k = 3
    args.gmt_topk_list = "3"
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.emit_small_gene_sets = True
    drug_response_screen.run(args)
    captured = capsys.readouterr()
    assert "broad receptor-family targets" in captured.err


def test_drug_response_bundle_resources_are_used(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    response = tmp_path / "resp.tsv"
    response.write_text(
        "\n".join(
            [
                "sample_id\tdrug_id\tresponse",
                "S1\tDrug A\t0.8",
                "S1\tDrug B\t0.8",
                "S1\tDrug Warn\t0.9",
                "S1\tDrug Drop\t0.7",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    args = Args()
    args.out_dir = str(tmp_path / "bundle_run")
    args.response_tsv = str(response)
    args.drug_targets_tsv = None
    args.groups_tsv = None
    args.contrast_method = "none"
    args.response_metric = "other"
    args.response_direction = "higher_is_more_sensitive"
    args.response_transform = "none"
    args.max_programs = 1
    args.top_k = 2
    args.gmt_topk_list = "2"
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.emit_small_gene_sets = True
    args.resources_dir = "tests/data/drug_response_bundle"
    result = drug_response_screen.run(args)
    assert result["out_dir"] == str(out_dir := Path(args.out_dir))
    row = _read_manifest_rows(out_dir)[0]
    genes = _read_geneset_rows(out_dir / row["path"] / "geneset.tsv")
    assert genes[0]["gene_id"] == "GENESEL"
    meta = json.loads((out_dir / row["path"] / "geneset.meta.json").read_text(encoding="utf-8"))
    resources = meta["summary"]["target_summary"]["resources"]
    used_ids = {entry["id"] for entry in resources["used"]}
    assert "drug_target_edges_human_v1" in used_ids
    assert "drug_alias_map_human_v1" in used_ids
    assert "target_ubiquity_human_v1" in used_ids
    assert "compound_qc_human_v1" in used_ids
    root_summary = json.loads((out_dir / "run_summary.json").read_text(encoding="utf-8"))
    assert root_summary["target_summary"]["target_source"] == "bundle:drug_target_edges_human_v1"
    captured = capsys.readouterr()
    assert "recommended_use=warn" in captured.err


def test_drug_response_prism_fallback_without_bundle_targets_warns(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    args = Args()
    args.out_dir = str(tmp_path / "prism_fallback")
    args.response_tsv = None
    args.drug_targets_tsv = None
    args.groups_tsv = None
    args.contrast_method = None
    args.prism_matrix_csv = "tests/data/toy_prism_matrix.csv"
    args.prism_treatment_info_csv = "tests/data/toy_prism_treatment_info.csv"
    args.prism_cell_line_info_csv = "tests/data/toy_prism_cell_line_info.csv"
    args.resources_dir = str(tmp_path / "empty_resources")
    Path(args.resources_dir).mkdir(parents=True, exist_ok=True)
    result = drug_response_screen.run(args)
    assert result["n_groups"] == 2
    out_dir = Path(args.out_dir)
    summary = json.loads((out_dir / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["target_summary"]["target_source"] == "prism_treatment_info"
    captured = capsys.readouterr()
    assert "bundle_resource_missing resource=drug_target_edges_human_v1" in captured.err


def test_drug_response_metadata_records_missing_bundle_resources(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "missing_bundle")
    args.groups_tsv = None
    args.contrast_method = "none"
    args.resources_dir = str(tmp_path / "empty_resources")
    Path(args.resources_dir).mkdir(parents=True, exist_ok=True)
    drug_response_screen.run(args)
    out_dir = Path(args.out_dir)
    row = _read_manifest_rows(out_dir)[0]
    meta = json.loads((out_dir / row["path"] / "geneset.meta.json").read_text(encoding="utf-8"))
    resources = meta["summary"]["target_summary"]["resources"]
    missing_ids = {entry["resource_id"] for entry in resources["missing"]}
    assert "drug_alias_map_human_v1" in missing_ids
    assert "target_ubiquity_human_v1" in missing_ids
