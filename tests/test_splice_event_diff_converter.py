import csv
import gzip
import json
import shutil
from pathlib import Path

from geneset_extractors.converters import splice_event_diff
from geneset_extractors.core.validate import validate_output_dir


class Args:
    splice_tsv = "tests/data/toy_splice_event_diff.tsv"
    cluster_stats_tsv = None
    out_dir = "tests/tmp/splice_event_diff"
    organism = "human"
    genome_build = "human"
    signature_name = "toy_splice"
    dataset_label = "toy_splice_dataset"
    tool_family = "generic"
    event_id_column = None
    event_group_column = None
    event_type_column = None
    gene_id_column = None
    gene_symbol_column = None
    chrom_column = None
    start_column = None
    end_column = None
    strand_column = None
    score_column = None
    stat_column = None
    delta_psi_column = None
    psi_column = None
    padj_column = None
    pvalue_column = None
    probability_column = None
    read_support_column = None
    novel_flag_column = None
    annotation_status_column = None
    score_mode = "auto"
    score_transform = "signed"
    confidence_weight_mode = "combined"
    min_probability = 0.8
    min_read_support = 10.0
    neglog10p_cap = 50.0
    neglog10p_eps = 1e-300
    delta_psi_soft_floor = 0.05
    delta_psi_soft_floor_mode = "auto"
    event_dup_policy = "highest_confidence"
    gene_aggregation = "signed_topk_mean"
    gene_topk_events = 3
    gene_burden_penalty_mode = "auto"
    min_gene_burden_penalty = 0.35
    ambiguous_gene_policy = "drop"
    impact_mode = "conservative"
    impact_min = 0.75
    impact_max = 1.35
    resources_manifest = None
    resources_dir = None
    resource_policy = "skip"
    use_reference_bundle = True
    event_alias_resource_id = None
    event_ubiquity_resource_id = None
    event_impact_resource_id = None
    gene_burden_resource_id = None
    select = "top_k"
    top_k = 200
    quantile = 0.01
    min_score = 0.0
    normalize = "within_set_l1"
    emit_full = True
    emit_gmt = True
    gmt_out = None
    gmt_prefer_symbol = True
    gmt_require_symbol = True
    gmt_biotype_allowlist = "protein_coding"
    gmt_min_genes = 1
    gmt_max_genes = 10
    gmt_topk_list = "2"
    gmt_mass_list = ""
    gmt_split_signed = True
    gmt_format = "dig2col"
    emit_small_gene_sets = True


def _copy_resources(resources_dir: Path) -> None:
    resources_dir.mkdir(parents=True, exist_ok=True)
    for filename in (
        "splice_event_aliases_human_v1.tsv.gz",
        "splice_event_ubiquity_human_v1.tsv.gz",
        "splice_event_impact_human_v1.tsv.gz",
        "splice_gene_event_burden_human_v1.tsv.gz",
    ):
        shutil.copyfile(Path("tests/data") / filename, resources_dir / filename)


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _score_map(path: Path) -> dict[str, float]:
    return {row["gene_id"]: float(row["score"]) for row in _read_rows(path)}


def test_splice_event_diff_end_to_end_with_bundle(tmp_path: Path):
    resources_dir = tmp_path / "resources"
    _copy_resources(resources_dir)

    args = Args()
    args.out_dir = str(tmp_path / "splice_run")
    args.resources_dir = str(resources_dir)
    result = splice_event_diff.run(args)
    assert result["resolved_score_mode"] == "delta_psi_times_neglog10p"

    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)

    rows = _read_rows(Path(args.out_dir) / "geneset.tsv")
    assert rows
    assert rows[0]["gene_symbol"] == "KCNN4"
    assert "gene_burden_penalty" in rows[0]
    assert abs(sum(float(r["weight"]) for r in rows) - 1.0) < 1e-9

    gmt_lines = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8").strip().splitlines()
    assert any("__pos__" in line for line in gmt_lines)
    assert any("__neg__" in line for line in gmt_lines)

    meta = json.loads((Path(args.out_dir) / "geneset.meta.json").read_text(encoding="utf-8"))
    resources_info = meta["converter"]["parameters"]["resources"]
    assert {item["id"] for item in resources_info["used"]} == {
        "splice_event_aliases_human_v1",
        "splice_event_ubiquity_human_v1",
        "splice_event_impact_human_v1",
        "splice_gene_event_burden_human_v1",
    }


def test_splice_event_diff_runs_without_bundle(tmp_path: Path, capsys):
    args = Args()
    args.out_dir = str(tmp_path / "splice_no_bundle")
    args.resources_dir = None
    args.use_reference_bundle = False
    result = splice_event_diff.run(args)
    assert result["n_genes"] >= 1
    captured = capsys.readouterr()
    assert "Missing" not in captured.err


def test_splice_event_diff_impact_changes_scores(tmp_path: Path):
    resources_dir = tmp_path / "resources"
    _copy_resources(resources_dir)

    args_imp = Args()
    args_imp.out_dir = str(tmp_path / "splice_imp")
    args_imp.resources_dir = str(resources_dir)
    args_imp.impact_mode = "conservative"
    splice_event_diff.run(args_imp)

    args_none = Args()
    args_none.out_dir = str(tmp_path / "splice_none")
    args_none.resources_dir = str(resources_dir)
    args_none.impact_mode = "none"
    splice_event_diff.run(args_none)

    imp_scores = _score_map(Path(args_imp.out_dir) / "geneset.full.tsv")
    none_scores = _score_map(Path(args_none.out_dir) / "geneset.full.tsv")
    assert imp_scores["G_KCNN4"] != none_scores["G_KCNN4"]


def test_splice_event_diff_bundle_matched_fixture_changes_scores(tmp_path: Path):
    resources_dir = tmp_path / "resources"
    _copy_resources(resources_dir)

    args_bundle = Args()
    args_bundle.splice_tsv = "tests/data/toy_splice_event_diff_matched.tsv"
    args_bundle.out_dir = str(tmp_path / "matched_bundle")
    args_bundle.resources_dir = str(resources_dir)
    splice_event_diff.run(args_bundle)

    args_none = Args()
    args_none.splice_tsv = "tests/data/toy_splice_event_diff_matched.tsv"
    args_none.out_dir = str(tmp_path / "matched_none")
    args_none.use_reference_bundle = False
    splice_event_diff.run(args_none)

    bundle_scores = _score_map(Path(args_bundle.out_dir) / "geneset.full.tsv")
    none_scores = _score_map(Path(args_none.out_dir) / "geneset.full.tsv")
    assert bundle_scores != none_scores
    assert bundle_scores["G_MAPK1"] != none_scores["G_MAPK1"]
    summary = json.loads((Path(args_bundle.out_dir) / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["n_events_with_low_confidence_prior_match"] >= 1


def test_splice_event_diff_low_confidence_bundle_priors_are_neutralized(tmp_path: Path):
    resources_dir = tmp_path / "resources"
    _copy_resources(resources_dir)

    args_bundle = Args()
    args_bundle.splice_tsv = "tests/data/toy_splice_event_diff_matched.tsv"
    args_bundle.out_dir = str(tmp_path / "matched_bundle_neutral")
    args_bundle.resources_dir = str(resources_dir)
    args_bundle.gene_burden_penalty_mode = "none"
    splice_event_diff.run(args_bundle)

    args_none = Args()
    args_none.splice_tsv = "tests/data/toy_splice_event_diff_matched.tsv"
    args_none.out_dir = str(tmp_path / "matched_none_neutral")
    args_none.use_reference_bundle = False
    args_none.gene_burden_penalty_mode = "none"
    splice_event_diff.run(args_none)

    bundle_scores = _score_map(Path(args_bundle.out_dir) / "geneset.full.tsv")
    none_scores = _score_map(Path(args_none.out_dir) / "geneset.full.tsv")
    assert abs(bundle_scores["G_MAPK1"] - none_scores["G_MAPK1"]) < 0.002
    assert bundle_scores["G_KCNN4"] != none_scores["G_KCNN4"]
    summary = json.loads((Path(args_bundle.out_dir) / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["n_low_confidence_ubiquity_priors_neutralized"] >= 1
    assert summary["n_low_confidence_impact_priors_neutralized"] >= 1


def test_splice_event_diff_gene_burden_penalty_changes_scores(tmp_path: Path):
    resources_dir = tmp_path / "resources"
    _copy_resources(resources_dir)
    with gzip.open(resources_dir / "splice_gene_event_burden_human_v1.tsv.gz", "wt", encoding="utf-8", newline="") as fh:
        fh.write(
            "\n".join(
                [
                    "gene_symbol\tn_canonical_events_ref\tn_high_confidence_events_ref\tn_low_confidence_events_ref\tn_unique_event_groups_ref\tn_studies_ref\tn_studies_high_confidence_ref\tfraction_low_confidence_events_ref\tmedian_unique_groups_per_study",
                    "GENEA\t8\t8\t0\t6\t4\t4\t0.0\t2.0",
                    "GENEB\t1\t1\t0\t1\t1\t1\t0.0\t1.0",
                ]
            )
            + "\n"
        )

    args_none = Args()
    args_none.splice_tsv = "tests/data/toy_splice_burden_diff.tsv"
    args_none.out_dir = str(tmp_path / "burden_none")
    args_none.use_reference_bundle = False
    args_none.gene_aggregation = "sum"
    args_none.gene_topk_events = 1
    args_none.gene_burden_penalty_mode = "none"
    splice_event_diff.run(args_none)

    args_pen = Args()
    args_pen.splice_tsv = "tests/data/toy_splice_burden_diff.tsv"
    args_pen.out_dir = str(tmp_path / "burden_current")
    args_pen.resources_dir = str(resources_dir)
    args_pen.gene_aggregation = "sum"
    args_pen.gene_topk_events = 1
    args_pen.gene_burden_penalty_mode = "reference_bundle"
    splice_event_diff.run(args_pen)

    none_scores = _score_map(Path(args_none.out_dir) / "geneset.full.tsv")
    pen_scores = _score_map(Path(args_pen.out_dir) / "geneset.full.tsv")
    assert abs(pen_scores["G_A"]) < abs(none_scores["G_A"])
    pen_rows = _read_rows(Path(args_pen.out_dir) / "geneset.tsv")
    assert pen_rows[0]["gene_symbol"] == "GENEB"


def test_splice_event_diff_event_group_collapse_reduces_cluster_stacking(tmp_path: Path):
    args = Args()
    args.splice_tsv = "tests/data/toy_splice_burden_diff.tsv"
    args.out_dir = str(tmp_path / "event_group_collapse")
    args.use_reference_bundle = False
    args.gene_aggregation = "sum"
    args.gene_burden_penalty_mode = "none"
    splice_event_diff.run(args)

    rows = _read_rows(Path(args.out_dir) / "geneset.tsv")
    assert rows[0]["gene_symbol"] == "GENEB"
    summary = json.loads((Path(args.out_dir) / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["n_independent_event_groups"] == 2


def test_splice_event_diff_ambiguous_gene_handling(tmp_path: Path):
    args_drop = Args()
    args_drop.out_dir = str(tmp_path / "splice_drop")
    args_drop.use_reference_bundle = False
    args_drop.ambiguous_gene_policy = "drop"
    splice_event_diff.run(args_drop)
    genes_drop = {row["gene_symbol"] for row in _read_rows(Path(args_drop.out_dir) / "geneset.full.tsv")}
    assert "GENEA" not in genes_drop
    assert "GENEB" not in genes_drop

    args_split = Args()
    args_split.out_dir = str(tmp_path / "splice_split")
    args_split.use_reference_bundle = False
    args_split.ambiguous_gene_policy = "split_equal"
    splice_event_diff.run(args_split)
    genes_split = {row["gene_symbol"] for row in _read_rows(Path(args_split.out_dir) / "geneset.full.tsv")}
    assert {"GENEA", "GENEB"}.issubset(genes_split)


def test_splice_event_diff_duplicate_event_collapse(tmp_path: Path):
    args_conf = Args()
    args_conf.out_dir = str(tmp_path / "splice_conf")
    args_conf.use_reference_bundle = False
    args_conf.event_dup_policy = "highest_confidence"
    splice_event_diff.run(args_conf)

    args_abs = Args()
    args_abs.out_dir = str(tmp_path / "splice_abs")
    args_abs.use_reference_bundle = False
    args_abs.event_dup_policy = "max_abs"
    splice_event_diff.run(args_abs)

    conf_scores = _score_map(Path(args_conf.out_dir) / "geneset.full.tsv")
    abs_scores = _score_map(Path(args_abs.out_dir) / "geneset.full.tsv")
    assert conf_scores["G_KCNN4"] != abs_scores["G_KCNN4"] or conf_scores["G_KCNN4"] > 0


def test_splice_event_diff_warning_path(tmp_path: Path, capsys):
    warning_path = tmp_path / "warning.tsv"
    warning_path.write_text(
        "\n".join(
            [
                "event_id\tevent_group\tevent_type\tgene_id\tgene_symbol\tdelta_psi\tpvalue\tprobability\tread_support\tannotation_status",
                "W1\tG1\tretained_intron\tG_R1\tR1\t0.04\t0.2\t0.4\t1\tnovel",
                "W2\tG2\tretained_intron\tG_R2\tR2\t0.03\t0.3\t0.3\t1\tnovel",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    args = Args()
    args.splice_tsv = str(warning_path)
    args.out_dir = str(tmp_path / "warning_run")
    args.use_reference_bundle = False
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "2"
    splice_event_diff.run(args)

    captured = capsys.readouterr()
    assert "retained introns dominate" in captured.err or "low-confidence events dominate" in captured.err
    summary = json.loads((Path(args.out_dir) / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["warnings"]
    assert summary["n_delta_psi_soft_floor_downweighted"] >= 1


def test_splice_event_diff_locus_density_warning_and_penalty(tmp_path: Path):
    locus_path = tmp_path / "locus.tsv"
    locus_path.write_text(
        "\n".join(
            [
                "event_id\tevent_group\tevent_type\tgene_id\tgene_symbol\tchrom\tstart\tend\tstrand\tdelta_psi\tpvalue\tprobability\tread_support\tannotation_status",
                "L1\tLG1\texon_skip\tG1\tGENE1\tchr1\t100000\t100120\t+\t0.60\t1e-6\t0.99\t50\tannotated_coding",
                "L2\tLG2\texon_skip\tG2\tGENE2\tchr1\t150000\t150120\t+\t0.58\t1e-6\t0.99\t50\tannotated_coding",
                "L3\tLG3\texon_skip\tG3\tGENE3\tchr1\t250000\t250120\t+\t0.56\t1e-6\t0.99\t50\tannotated_coding",
                "L4\tLG4\texon_skip\tG4\tGENE4\tchr1\t350000\t350120\t+\t0.54\t1e-6\t0.99\t50\tannotated_coding",
                "L5\tLG5\texon_skip\tG5\tGENE5\tchr1\t450000\t450120\t+\t0.52\t1e-6\t0.99\t50\tannotated_coding",
                "L6\tLG6\texon_skip\tG6\tGENE6\tchr8\t100000\t100120\t+\t0.40\t1e-6\t0.99\t50\tannotated_coding",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    args_none = Args()
    args_none.splice_tsv = str(locus_path)
    args_none.out_dir = str(tmp_path / "locus_none")
    args_none.use_reference_bundle = False
    args_none.top_k = 10
    args_none.locus_density_top_n = 5
    args_none.locus_density_penalty_mode = "none"
    splice_event_diff.run(args_none)

    args_pen = Args()
    args_pen.splice_tsv = str(locus_path)
    args_pen.out_dir = str(tmp_path / "locus_pen")
    args_pen.use_reference_bundle = False
    args_pen.top_k = 10
    args_pen.locus_density_top_n = 5
    args_pen.locus_density_penalty_mode = "window_diversity"
    splice_event_diff.run(args_pen)

    summary = json.loads((Path(args_none.out_dir) / "run_summary.json").read_text(encoding="utf-8"))
    assert "likely_locus_density_artifact" in summary["warnings"]

    none_scores = _score_map(Path(args_none.out_dir) / "geneset.full.tsv")
    pen_scores = _score_map(Path(args_pen.out_dir) / "geneset.full.tsv")
    assert pen_scores["G2"] != none_scores["G2"]
