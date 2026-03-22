import csv
import gzip
import json
import shutil
from pathlib import Path

from geneset_extractors.converters import splice_event_diff
from geneset_extractors.core.validate import validate_output_dir
from tests.provenance_helpers import (
    assert_node_has_structured_resource_metadata,
    file_node_for_role,
    load_provenance,
)


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
    gene_support_penalty_mode = "auto"
    locus_density_penalty_mode = "none"
    locus_density_penalty_mode_explicit = False
    locus_density_window_bp = 20000000
    locus_density_top_n = 20
    ambiguous_gene_policy = "drop"
    impact_mode = "conservative"
    impact_min = 0.75
    impact_max = 1.35
    bundle_prior_profile = "auto"
    bundle_prior_profile_explicit = False
    artifact_action = "warn"
    artifact_action_explicit = False
    resources_manifest = None
    resources_dir = None
    resource_policy = "skip"
    use_reference_bundle = True
    event_alias_resource_id = None
    event_ubiquity_resource_id = None
    event_ubiquity_by_dataset_resource_id = None
    event_impact_resource_id = None
    gene_burden_resource_id = None
    gene_burden_by_dataset_resource_id = None
    gene_locus_resource_id = None
    source_dataset = None
    bundle_same_dataset_policy = "exclude"
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
    provenance_overlay_json = None


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


def _summary(out_dir: str | Path) -> dict[str, object]:
    return json.loads((Path(out_dir) / "run_summary.json").read_text(encoding="utf-8"))


def _write_gz_tsv(path: Path, text: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        fh.write(text)


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
    provenance = load_provenance(args.out_dir)
    node = file_node_for_role(provenance, "event_alias_table")
    assert_node_has_structured_resource_metadata(node)


def test_splice_event_diff_runs_without_bundle(tmp_path: Path, capsys):
    args = Args()
    args.out_dir = str(tmp_path / "splice_no_bundle")
    args.resources_dir = None
    args.use_reference_bundle = False
    result = splice_event_diff.run(args)
    assert result["n_genes"] >= 1
    captured = capsys.readouterr()
    assert "Missing" not in captured.err


def test_splice_event_diff_tcga_dynamic_defaults_activate_only_when_needed(tmp_path: Path):
    dynamic_path = tmp_path / "dynamic.tsv"
    dynamic_path.write_text(
        "\n".join(
            [
                "event_id\tevent_group\tevent_type\tgene_id\tgene_symbol\tchrom\tstart\tend\tstrand\tdelta_psi\tpvalue\tprobability\tread_support\tannotation_status",
                "D1\tDG1\texon_skip\tG_KCNN4\tKCNN4\tchr19\t100\t180\t+\t0.50\t1e-6\t0.99\t30\tannotated_coding",
                "D2\tDG2\talt_donor\tG_MAPK1\tMAPK1\tchr22\t300\t390\t-\t0.40\t1e-5\t0.98\t25\tannotated_coding",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    args_generic = Args()
    args_generic.splice_tsv = str(dynamic_path)
    args_generic.out_dir = str(tmp_path / "generic_dynamic")
    args_generic.use_reference_bundle = False
    args_generic.tool_family = "generic"
    splice_event_diff.run(args_generic)
    generic_summary = _summary(args_generic.out_dir)
    assert generic_summary["tcga_public_mode"] is False
    assert generic_summary["locus_density_penalty_mode_effective"] == "none"
    assert generic_summary["effective_artifact_action"] == "warn"
    assert generic_summary["effective_bundle_prior_profile"] == "none"

    args_tcga = Args()
    args_tcga.splice_tsv = str(dynamic_path)
    args_tcga.out_dir = str(tmp_path / "tcga_dynamic")
    args_tcga.use_reference_bundle = False
    args_tcga.tool_family = "tcga_spliceseq"
    splice_event_diff.run(args_tcga)
    tcga_summary = _summary(args_tcga.out_dir)
    assert tcga_summary["tcga_public_mode"] is True
    assert tcga_summary["locus_density_penalty_mode_effective"] == "auto"
    assert tcga_summary["effective_artifact_action"] == "suppress_if_persistent"
    assert tcga_summary["effective_bundle_prior_profile"] == "none"

    args_override = Args()
    args_override.splice_tsv = str(dynamic_path)
    args_override.out_dir = str(tmp_path / "tcga_override")
    args_override.use_reference_bundle = False
    args_override.tool_family = "tcga_spliceseq"
    args_override.locus_density_penalty_mode = "chromosome_diversity"
    args_override.locus_density_penalty_mode_explicit = True
    args_override.artifact_action = "prune"
    args_override.artifact_action_explicit = True
    args_override.bundle_prior_profile = "full"
    args_override.bundle_prior_profile_explicit = True
    splice_event_diff.run(args_override)
    override_summary = _summary(args_override.out_dir)
    assert override_summary["locus_density_penalty_mode_effective"] == "chromosome_diversity"
    assert override_summary["effective_artifact_action"] == "prune"
    assert override_summary["effective_bundle_prior_profile"] == "full"


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


def test_splice_event_diff_auto_prior_profile_degrades_medium_self_only_but_keeps_high_confidence_full(tmp_path: Path):
    resources_dir = tmp_path / "resources_auto"
    resources_dir.mkdir(parents=True, exist_ok=True)
    _write_gz_tsv(
        resources_dir / "splice_event_aliases_human_v1.tsv.gz",
        "\n".join(
            [
                "input_event_key\tcanonical_event_key\tcanonicalization_status\tcanonicalization_confidence\tevent_key_namespace\tgene_id\tgene_symbol\tevent_type\tchrom\tstart\tend\tstrand\tsource_dataset",
                "M1\ttcga_spliceseq_asid::KCNN4::exon_skip::M1\tsource_family_stable_id\tmedium\ttcga_spliceseq_asid\tG_KCNN4\tKCNN4\texon_skip\t\t\t\t\tTCGA_BRCA",
            ]
        )
        + "\n",
    )
    _write_gz_tsv(
        resources_dir / "splice_event_ubiquity_human_v1.tsv.gz",
        "\n".join(
            [
                "canonical_event_key\tgene_id\tgene_symbol\tevent_type\tcanonicalization_status\tcanonicalization_confidence\tevent_key_namespace\tn_samples_ref\tdf_ref\tfraction_ref\tidf_ref\tn_datasets_ref\tfraction_datasets_ref\tn_source_datasets_ref\tmax_dataset_fraction_ref",
                "tcga_spliceseq_asid::KCNN4::exon_skip::M1\tG_KCNN4\tKCNN4\texon_skip\tsource_family_stable_id\tmedium\ttcga_spliceseq_asid\t5\t5\t1.0\t1.0\t1\t1.0\t1\t1.0",
                "chr7:300-360:+:alt_donor:HG1\tG_MAPK1\tMAPK1\talt_donor\tcoordinate_canonical\thigh\tglobal_coordinate\t12\t4\t0.33\t1.8\t3\t1.0\t3\t0.5",
            ]
        )
        + "\n",
    )
    _write_gz_tsv(
        resources_dir / "splice_event_ubiquity_by_dataset_human_v1.tsv.gz",
        "\n".join(
            [
                "canonical_event_key\tsource_dataset\tgene_id\tgene_symbol\tevent_type\tcanonicalization_status\tcanonicalization_confidence\tevent_key_namespace\tn_samples_ref\tdf_ref",
                "tcga_spliceseq_asid::KCNN4::exon_skip::M1\tTCGA_BRCA\tG_KCNN4\tKCNN4\texon_skip\tsource_family_stable_id\tmedium\ttcga_spliceseq_asid\t5\t5",
                "chr7:300-360:+:alt_donor:HG1\tTCGA_LUAD\tG_MAPK1\tMAPK1\talt_donor\tcoordinate_canonical\thigh\tglobal_coordinate\t4\t2",
                "chr7:300-360:+:alt_donor:HG1\tTCGA_HNSC\tG_MAPK1\tMAPK1\talt_donor\tcoordinate_canonical\thigh\tglobal_coordinate\t4\t1",
            ]
        )
        + "\n",
    )
    _write_gz_tsv(
        resources_dir / "splice_event_impact_human_v1.tsv.gz",
        "\n".join(
            [
                "canonical_event_key\tgene_id\tgene_symbol\tevent_type\tcanonicalization_status\tcanonicalization_confidence\tevent_key_namespace\timpact_weight_raw\timpact_evidence\tannotation_status\tn_datasets_ref\tn_source_datasets_ref\tmax_dataset_fraction_ref",
                "tcga_spliceseq_asid::KCNN4::exon_skip::M1\tG_KCNN4\tKCNN4\texon_skip\tsource_family_stable_id\tmedium\ttcga_spliceseq_asid\t1.2\t0.6\tannotated_coding\t1\t1\t1.0",
                "chr7:300-360:+:alt_donor:HG1\tG_MAPK1\tMAPK1\talt_donor\tcoordinate_canonical\thigh\tglobal_coordinate\t1.2\t0.6\tannotated_coding\t3\t3\t0.5",
            ]
        )
        + "\n",
    )
    _write_gz_tsv(
        resources_dir / "splice_gene_event_burden_human_v1.tsv.gz",
        "\n".join(
            [
                "gene_symbol\tn_canonical_events_ref\tn_high_confidence_events_ref\tn_medium_confidence_events_ref\tn_low_confidence_events_ref\tn_unique_event_groups_ref\tn_studies_ref\tn_studies_high_confidence_ref\tfraction_low_confidence_events_ref\tmedian_unique_groups_per_study\tmax_dataset_fraction_ref",
                "KCNN4\t1\t0\t1\t0\t1\t1\t0\t0.0\t1.0\t1.0",
                "MAPK1\t1\t1\t0\t0\t1\t3\t3\t0.0\t1.0\t0.5",
            ]
        )
        + "\n",
    )
    _write_gz_tsv(
        resources_dir / "splice_gene_event_burden_by_dataset_human_v1.tsv.gz",
        "\n".join(
            [
                "gene_symbol\tsource_dataset\tn_canonical_events_ref\tn_high_confidence_events_ref\tn_medium_confidence_events_ref\tn_low_confidence_events_ref\tn_unique_event_groups_ref",
                "KCNN4\tTCGA_BRCA\t1\t0\t1\t0\t1",
                "MAPK1\tTCGA_LUAD\t1\t1\t0\t0\t1",
                "MAPK1\tTCGA_HNSC\t1\t1\t0\t0\t1",
            ]
        )
        + "\n",
    )

    medium_path = tmp_path / "medium_tcga.tsv"
    medium_path.write_text(
        "\n".join(
            [
                "event_id\tevent_group\tevent_type\tgene_id\tgene_symbol\tdelta_psi\tpvalue\tprobability\tread_support\tannotation_status",
                "M1\tMG1\texon_skip\tG_KCNN4\tKCNN4\t0.75\t1e-6\t0.99\t30\tannotated_coding",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    args_medium = Args()
    args_medium.splice_tsv = str(medium_path)
    args_medium.out_dir = str(tmp_path / "medium_runtime")
    args_medium.tool_family = "tcga_spliceseq"
    args_medium.resources_dir = str(resources_dir)
    args_medium.source_dataset = "TCGA_BRCA"
    splice_event_diff.run(args_medium)
    medium_summary = _summary(args_medium.out_dir)
    assert medium_summary["effective_bundle_prior_profile"] == "none"
    assert medium_summary["bundle_prior_downgrade_reason"] == "no_meaningful_nonself_bundle_support"
    assert medium_summary["n_same_dataset_ubiquity_priors_neutralized"] >= 1

    high_path = tmp_path / "high_tcga.tsv"
    high_path.write_text(
        "\n".join(
            [
                "event_id\tevent_group\tevent_type\tgene_id\tgene_symbol\tchrom\tstart\tend\tstrand\tdelta_psi\tpvalue\tprobability\tread_support\tannotation_status",
                "H1\tHG1\talt_donor\tG_MAPK1\tMAPK1\tchr7\t300\t360\t+\t0.65\t1e-6\t0.99\t30\tannotated_coding",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    args_high = Args()
    args_high.splice_tsv = str(high_path)
    args_high.out_dir = str(tmp_path / "high_runtime")
    args_high.tool_family = "tcga_spliceseq"
    args_high.resources_dir = str(resources_dir)
    args_high.source_dataset = "TCGA_BRCA"
    splice_event_diff.run(args_high)
    high_summary = _summary(args_high.out_dir)
    assert high_summary["effective_bundle_prior_profile"] == "full"
    assert high_summary["bundle_prior_downgrade_reason"] == "sufficient_nonself_support"
    assert high_summary["n_events_matched_to_ubiquity_prior"] == 1


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


def test_splice_event_diff_tcga_persistent_locus_artifact_is_withheld(tmp_path: Path):
    locus_path = tmp_path / "tcga_artifact.tsv"
    locus_path.write_text(
        "\n".join(
            [
                "event_id\tevent_group\tevent_type\tgene_id\tgene_symbol\tchrom\tstart\tend\tstrand\tdelta_psi\tpvalue\tprobability\tread_support\tannotation_status",
                "T1\tTG1\texon_skip\tG1\tGENE1\tchr1\t100000\t100120\t+\t0.70\t1e-8\t0.99\t50\tannotated_coding",
                "T2\tTG2\texon_skip\tG2\tGENE2\tchr1\t200000\t200120\t+\t0.69\t1e-8\t0.99\t50\tannotated_coding",
                "T3\tTG3\texon_skip\tG3\tGENE3\tchr1\t300000\t300120\t+\t0.68\t1e-8\t0.99\t50\tannotated_coding",
                "T4\tTG4\texon_skip\tG4\tGENE4\tchr1\t400000\t400120\t+\t0.67\t1e-8\t0.99\t50\tannotated_coding",
                "T5\tTG5\texon_skip\tG5\tGENE5\tchr1\t500000\t500120\t+\t0.66\t1e-8\t0.99\t50\tannotated_coding",
                "T6\tTG6\texon_skip\tG6\tGENE6\tchr1\t600000\t600120\t+\t0.65\t1e-8\t0.99\t50\tannotated_coding",
                "T7\tTG7\texon_skip\tG7\tGENE7\tchr8\t100000\t100120\t+\t0.40\t1e-8\t0.99\t50\tannotated_coding",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    args = Args()
    args.splice_tsv = str(locus_path)
    args.out_dir = str(tmp_path / "tcga_artifact_run")
    args.use_reference_bundle = False
    args.tool_family = "tcga_spliceseq"
    args.top_k = 10
    args.locus_density_top_n = 5
    splice_event_diff.run(args)

    summary = _summary(args.out_dir)
    assert summary["tcga_public_mode"] is True
    assert summary["locus_density_penalty_mode_effective"] == "auto"
    assert summary["effective_artifact_action"] == "suppress_if_persistent"
    assert summary["quality_tier"] == "artifact_prone_withheld"
    assert summary["artifact_action_applied"] == "withheld"
    assert summary["output_withheld_reason"] == "persistent_locus_density_artifact"
    assert "persistent_locus_density_artifact" in summary["artifact_flags"]
    assert summary["artifact_summary"]["likely_persistent_artifact"] is True
    assert summary["locus_density_after"]["likely_artifact"] is True
    assert _read_rows(Path(args.out_dir) / "geneset.tsv") == []
    assert _read_rows(Path(args.out_dir) / "geneset.full.tsv")
    assert not (Path(args.out_dir) / "genesets.gmt").exists()


def test_splice_event_diff_generic_raw_id_is_not_promoted_to_medium(tmp_path: Path):
    raw_path = tmp_path / "generic_raw.tsv"
    raw_path.write_text(
        "\n".join(
            [
                "event_id\tevent_group\tevent_type\tgene_id\tgene_symbol\tdelta_psi\tpvalue\tprobability\tread_support\tannotation_status",
                "A123\tCG1\texon_skip\tG_KCNN4\tKCNN4\t0.8\t1e-6\t0.99\t30\tannotated_coding",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    args = Args()
    args.splice_tsv = str(raw_path)
    args.out_dir = str(tmp_path / "generic_raw_run")
    args.use_reference_bundle = False
    splice_event_diff.run(args)
    summary = json.loads((Path(args.out_dir) / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["canonicalization_confidence_counts"].get("medium", 0) == 0


def test_splice_event_diff_medium_priors_are_namespace_safe(tmp_path: Path):
    resources_dir = tmp_path / "resources_ns"
    resources_dir.mkdir(parents=True, exist_ok=True)
    _write_gz_tsv(
        resources_dir / "splice_event_aliases_human_v1.tsv.gz",
        "\n".join(
            [
                "input_event_key\tcanonical_event_key\tcanonicalization_status\tcanonicalization_confidence\tevent_key_namespace\tgene_id\tgene_symbol\tevent_type\tchrom\tstart\tend\tstrand\tsource_dataset",
                "A123\ttcga_spliceseq_asid::KCNN4::exon_skip::A123\tsource_family_stable_id\tmedium\ttcga_spliceseq_asid\tG_KCNN4\tKCNN4\texon_skip\t\t\t\t\tTCGA_X",
            ]
        )
        + "\n",
    )
    _write_gz_tsv(
        resources_dir / "splice_event_ubiquity_human_v1.tsv.gz",
        "\n".join(
            [
                "canonical_event_key\tgene_id\tgene_symbol\tevent_type\tcanonicalization_status\tcanonicalization_confidence\tevent_key_namespace\tn_samples_ref\tdf_ref\tfraction_ref\tidf_ref\tn_datasets_ref\tfraction_datasets_ref",
                "tcga_spliceseq_asid::KCNN4::exon_skip::A123\tG_KCNN4\tKCNN4\texon_skip\tsource_family_stable_id\tmedium\ttcga_spliceseq_asid\t10\t4\t0.4\t2.0\t2\t1.0",
            ]
        )
        + "\n",
    )
    _write_gz_tsv(
        resources_dir / "splice_event_ubiquity_by_dataset_human_v1.tsv.gz",
        "\n".join(
            [
                "canonical_event_key\tsource_dataset\tgene_id\tgene_symbol\tevent_type\tcanonicalization_status\tcanonicalization_confidence\tevent_key_namespace\tn_samples_ref\tdf_ref",
                "tcga_spliceseq_asid::KCNN4::exon_skip::A123\tTCGA_X\tG_KCNN4\tKCNN4\texon_skip\tsource_family_stable_id\tmedium\ttcga_spliceseq_asid\t5\t2",
                "tcga_spliceseq_asid::KCNN4::exon_skip::A123\tTCGA_Y\tG_KCNN4\tKCNN4\texon_skip\tsource_family_stable_id\tmedium\ttcga_spliceseq_asid\t5\t2",
            ]
        )
        + "\n",
    )
    _write_gz_tsv(
        resources_dir / "splice_event_impact_human_v1.tsv.gz",
        "\n".join(
            [
                "canonical_event_key\tgene_id\tgene_symbol\tevent_type\tcanonicalization_status\tcanonicalization_confidence\tevent_key_namespace\timpact_weight_raw\timpact_evidence\tannotation_status\tn_datasets_ref",
                "tcga_spliceseq_asid::KCNN4::exon_skip::A123\tG_KCNN4\tKCNN4\texon_skip\tsource_family_stable_id\tmedium\ttcga_spliceseq_asid\t1.2\t0.6\tannotated_coding\t2",
            ]
        )
        + "\n",
    )
    _write_gz_tsv(
        resources_dir / "splice_gene_event_burden_human_v1.tsv.gz",
        "gene_symbol\tn_canonical_events_ref\tn_high_confidence_events_ref\tn_medium_confidence_events_ref\tn_low_confidence_events_ref\tn_unique_event_groups_ref\tn_studies_ref\tn_studies_high_confidence_ref\tfraction_low_confidence_events_ref\tmedian_unique_groups_per_study\nKCNN4\t1\t0\t1\t0\t1\t2\t0\t0.0\t1.0\n",
    )
    _write_gz_tsv(
        resources_dir / "splice_gene_event_burden_by_dataset_human_v1.tsv.gz",
        "gene_symbol\tsource_dataset\tn_canonical_events_ref\tn_high_confidence_events_ref\tn_medium_confidence_events_ref\tn_low_confidence_events_ref\tn_unique_event_groups_ref\nKCNN4\tTCGA_X\t1\t0\t1\t0\t1\nKCNN4\tTCGA_Y\t1\t0\t1\t0\t1\n",
    )

    generic_path = tmp_path / "generic.tsv"
    generic_path.write_text(
        "\n".join(
            [
                "event_id\tevent_group\tevent_type\tgene_id\tgene_symbol\tdelta_psi\tpvalue\tprobability\tread_support\tannotation_status",
                "A123\tCG1\texon_skip\tG_KCNN4\tKCNN4\t0.8\t1e-6\t0.99\t30\tannotated_coding",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    tcga_path = tmp_path / "tcga.tsv"
    tcga_path.write_text(generic_path.read_text(encoding="utf-8"), encoding="utf-8")

    args_generic = Args()
    args_generic.splice_tsv = str(generic_path)
    args_generic.out_dir = str(tmp_path / "generic_bundle")
    args_generic.resources_dir = str(resources_dir)
    args_generic.tool_family = "generic"
    splice_event_diff.run(args_generic)

    args_tcga = Args()
    args_tcga.splice_tsv = str(tcga_path)
    args_tcga.out_dir = str(tmp_path / "tcga_bundle")
    args_tcga.resources_dir = str(resources_dir)
    args_tcga.tool_family = "tcga_spliceseq"
    splice_event_diff.run(args_tcga)

    generic_scores = _score_map(Path(args_generic.out_dir) / "geneset.full.tsv")
    tcga_scores = _score_map(Path(args_tcga.out_dir) / "geneset.full.tsv")
    assert tcga_scores["G_KCNN4"] != generic_scores["G_KCNN4"]
    generic_summary = json.loads((Path(args_generic.out_dir) / "run_summary.json").read_text(encoding="utf-8"))
    tcga_summary = json.loads((Path(args_tcga.out_dir) / "run_summary.json").read_text(encoding="utf-8"))
    assert generic_summary["canonicalization_confidence_counts"].get("medium", 0) == 0
    assert tcga_summary["canonicalization_confidence_counts"].get("medium", 0) == 1


def test_splice_event_diff_auto_locus_density_penalty_diversifies_top_genes(tmp_path: Path):
    locus_path = tmp_path / "locus_auto.tsv"
    locus_path.write_text(
        "\n".join(
            [
                "event_id\tevent_group\tevent_type\tgene_id\tgene_symbol\tchrom\tstart\tend\tstrand\tdelta_psi\tpvalue\tprobability\tread_support\tannotation_status",
                "L1\tLG1\texon_skip\tG1\tGENE1\tchr1\t100000\t100120\t+\t0.60\t1e-6\t0.99\t50\tannotated_coding",
                "L2\tLG2\texon_skip\tG2\tGENE2\tchr1\t150000\t150120\t+\t0.58\t1e-6\t0.99\t50\tannotated_coding",
                "L3\tLG3\texon_skip\tG3\tGENE3\tchr1\t250000\t250120\t+\t0.56\t1e-6\t0.99\t50\tannotated_coding",
                "L4\tLG4\texon_skip\tG4\tGENE4\tchr1\t350000\t350120\t+\t0.54\t1e-6\t0.99\t50\tannotated_coding",
                "L5\tLG5\texon_skip\tG5\tGENE5\tchr1\t450000\t450120\t+\t0.52\t1e-6\t0.99\t50\tannotated_coding",
                "L6\tLG6\texon_skip\tG6\tGENE6\tchr8\t100000\t100120\t+\t0.50\t1e-6\t0.99\t50\tannotated_coding",
                "L7\tLG7\texon_skip\tG7\tGENE7\tchr9\t100000\t100120\t+\t0.49\t1e-6\t0.99\t50\tannotated_coding",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    args = Args()
    args.splice_tsv = str(locus_path)
    args.out_dir = str(tmp_path / "locus_auto")
    args.use_reference_bundle = False
    args.top_k = 10
    args.locus_density_top_n = 5
    args.locus_density_penalty_mode = "auto"
    splice_event_diff.run(args)
    summary = json.loads((Path(args.out_dir) / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["locus_density"]["likely_artifact"] is True
    assert summary["locus_density"]["n_genes_penalized_for_locus_density"] >= 1
    rows = _read_rows(Path(args.out_dir) / "geneset.full.tsv")
    penalty_map = {row["gene_symbol"]: float(row["locus_density_penalty"]) for row in rows}
    assert penalty_map["GENE1"] == 1.0
    assert penalty_map["GENE2"] < 1.0


def test_splice_event_diff_tcga_nonartifact_signal_survives_dynamic_controls(tmp_path: Path):
    clean_path = tmp_path / "tcga_clean.tsv"
    clean_path.write_text(
        "\n".join(
            [
                "event_id\tevent_group\tevent_type\tgene_id\tgene_symbol\tchrom\tstart\tend\tstrand\tdelta_psi\tpvalue\tprobability\tread_support\tannotation_status",
                "C1\tCG1\texon_skip\tG1\tGENE1\tchr1\t100000\t100120\t+\t0.70\t1e-8\t0.99\t50\tannotated_coding",
                "C2\tCG2\talt_donor\tG2\tGENE2\tchr2\t200000\t200120\t+\t0.68\t1e-8\t0.99\t50\tannotated_coding",
                "C3\tCG3\texon_skip\tG3\tGENE3\tchr3\t300000\t300120\t+\t0.66\t1e-8\t0.99\t50\tannotated_coding",
                "C4\tCG4\talt_acceptor\tG4\tGENE4\tchr4\t400000\t400120\t+\t0.64\t1e-8\t0.99\t50\tannotated_coding",
                "C5\tCG5\texon_skip\tG5\tGENE5\tchr5\t500000\t500120\t+\t0.62\t1e-8\t0.99\t50\tannotated_coding",
                "C6\tCG6\talt_donor\tG6\tGENE6\tchr6\t600000\t600120\t+\t0.60\t1e-8\t0.99\t50\tannotated_coding",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    args = Args()
    args.splice_tsv = str(clean_path)
    args.out_dir = str(tmp_path / "tcga_clean_run")
    args.use_reference_bundle = False
    args.tool_family = "tcga_spliceseq"
    args.top_k = 10
    args.locus_density_top_n = 5
    splice_event_diff.run(args)
    summary = _summary(args.out_dir)
    assert summary["tcga_public_mode"] is True
    assert summary["locus_density_penalty_mode_effective"] == "auto"
    assert summary["effective_artifact_action"] == "suppress_if_persistent"
    assert summary["quality_tier"] == "standard"
    assert summary["output_withheld_reason"] == ""
    assert summary["artifact_summary"]["likely_persistent_artifact"] is False
    assert _read_rows(Path(args.out_dir) / "geneset.tsv")


def test_splice_event_diff_gene_support_penalty_rewards_coherent_multi_group_support(tmp_path: Path):
    support_path = tmp_path / "support.tsv"
    support_path.write_text(
        "\n".join(
            [
                "event_id\tevent_group\tevent_type\tgene_id\tgene_symbol\tdelta_psi\tpvalue\tprobability\tread_support\tannotation_status",
                "A1\tGA1\texon_skip\tG_A\tGENEA\t1.60\t1e-8\t0.99\t40\tannotated_coding",
                "A2\tGA2\texon_skip\tG_A\tGENEA\t-0.50\t1e-8\t0.99\t40\tannotated_coding",
                "B1\tGB1\texon_skip\tG_B\tGENEB\t0.55\t1e-8\t0.99\t40\tannotated_coding",
                "B2\tGB2\talt_donor\tG_B\tGENEB\t0.52\t1e-8\t0.99\t40\tannotated_coding",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    args_none = Args()
    args_none.splice_tsv = str(support_path)
    args_none.out_dir = str(tmp_path / "support_none")
    args_none.use_reference_bundle = False
    args_none.gene_support_penalty_mode = "none"
    args_none.gene_aggregation = "sum"
    splice_event_diff.run(args_none)

    args_auto = Args()
    args_auto.splice_tsv = str(support_path)
    args_auto.out_dir = str(tmp_path / "support_auto")
    args_auto.use_reference_bundle = False
    args_auto.gene_support_penalty_mode = "auto"
    args_auto.gene_aggregation = "sum"
    splice_event_diff.run(args_auto)

    none_rows = _read_rows(Path(args_none.out_dir) / "geneset.tsv")
    auto_rows = _read_rows(Path(args_auto.out_dir) / "geneset.tsv")
    assert none_rows[0]["gene_symbol"] == "GENEA"
    assert auto_rows[0]["gene_symbol"] == "GENEB"
    auto_full = _read_rows(Path(args_auto.out_dir) / "geneset.full.tsv")
    gene_b = next(row for row in auto_full if row["gene_symbol"] == "GENEB")
    assert int(gene_b["n_independent_event_groups_used"]) == 2
    assert float(gene_b["sign_coherence"]) > 0.9
