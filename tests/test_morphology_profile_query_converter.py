import csv
import json
from pathlib import Path

import pytest

from geneset_extractors.converters import morphology_profile_query
from geneset_extractors.core.validate import validate_output_dir


class Args:
    query_profiles_tsv = "tests/data/toy_morph_query_profiles.tsv"
    out_dir = "tests/tmp/morphology"
    organism = "human"
    genome_build = "hg38"
    dataset_label = "toy_morph"
    signature_name = "toy_query"

    query_id_column = "sample_id"
    query_profiles_delimiter = "\t"
    query_metadata_tsv = "tests/data/toy_morph_query_metadata.tsv"
    query_metadata_id_column = "sample_id"
    query_metadata_delimiter = "\t"
    group_query_by = "query_group"
    query_aggregate = "median"
    exclude_query_ids_from_reference = False

    reference_profiles_tsv = "tests/data/toy_morph_reference_profiles.tsv"
    reference_profiles_parquet = None
    reference_profiles_delimiter = "\t"
    reference_id_column = "perturbation_id"
    reference_metadata_tsv = "tests/data/toy_morph_reference_metadata.tsv"
    reference_metadata_id_column = "perturbation_id"
    reference_metadata_delimiter = "\t"
    compound_targets_tsv = "tests/data/toy_morph_compound_targets.tsv"
    compound_targets_delimiter = "\t"
    compound_id_column = "compound_id"
    compound_target_gene_symbol_column = "gene_symbol"
    compound_target_weight_column = "weight"
    feature_stats_tsv = "tests/data/toy_morph_feature_stats.tsv"
    feature_schema_tsv = "tests/data/toy_morph_feature_schema.tsv"

    resources_manifest = None
    resources_dir = None
    resource_policy = "skip"
    reference_bundle_id = None

    similarity_metric = "cosine"
    similarity_power = 1.0
    polarity = "both"
    max_reference_neighbors = 20
    min_similarity = 0.0
    hubness_penalty = "inverse_rank"
    min_specificity_confidence_to_emit_opposite = "medium"
    compound_weight = 0.5
    genetic_weight = 0.5

    select = "top_k"
    top_k = 2
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
    gmt_split_signed = False
    gmt_format = "classic"
    emit_small_gene_sets = True


def _manifest_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _geneset_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_morphology_profile_query_explicit_files(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "morph_explicit")
    result = morphology_profile_query.run(args)
    assert result["n_groups"] == 2

    out_dir = Path(args.out_dir)
    manifest = _manifest_rows(out_dir / "manifest.tsv")
    assert len(manifest) == 2
    assert (out_dir / "genesets.gmt").exists()

    q1_similar = out_dir / "program=Q1__polarity=similar" / "geneset.tsv"
    rows = _geneset_rows(q1_similar)
    assert [row["gene_id"] for row in rows] == ["GENE1", "GENE3"]
    weight_sum = sum(float(row["weight"]) for row in rows)
    assert abs(weight_sum - 1.0) < 1e-9

    meta = json.loads((out_dir / "program=Q1__polarity=similar" / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta["input"]["assay"] == "morphology"
    assert meta["summary"]["mapping_summary"]["n_compound_refs_with_targets"] == 2
    assert meta["gmt"]["written"] is True
    assert meta["summary"]["retrieval_confidence"] is not None
    assert meta["summary"]["specificity_confidence"] is not None
    assert isinstance(meta["summary"]["top_neighbor_ids"], list)

    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(out_dir, schema)


def test_morphology_profile_query_bundle_mode(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "morph_bundle")
    args.reference_profiles_tsv = None
    args.reference_metadata_tsv = None
    args.compound_targets_tsv = None
    args.feature_stats_tsv = None
    args.feature_schema_tsv = None
    args.resources_dir = "tests/data/morphology_bundle"
    args.reference_bundle_id = "morphology_jump_target_pilot_u2os_48h_v1"
    result = morphology_profile_query.run(args)
    assert result["n_groups"] == 2

    out_dir = Path(args.out_dir)
    meta = json.loads((out_dir / "program=Q1__polarity=similar" / "geneset.meta.json").read_text(encoding="utf-8"))
    resources = meta["summary"]["resources"]
    assert resources is not None
    used_ids = {row["id"] for row in resources["used"]}
    assert "morphology_jump_target_pilot_u2os_48h_v1" in used_ids
    assert meta["summary"]["parse_summary"]["bundle_manifest"]["_bundle_resolution"] == "local_resources_dir"


def test_morphology_profile_query_missing_bundle_warns_or_fails(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "morph_missing_bundle")
    args.reference_profiles_tsv = None
    args.reference_metadata_tsv = None
    args.compound_targets_tsv = None
    args.feature_stats_tsv = None
    args.feature_schema_tsv = None
    args.resources_dir = str(tmp_path / "empty_resources")
    args.reference_bundle_id = "morphology_jump_target_pilot_u2os_48h_v1"
    Path(args.resources_dir).mkdir(parents=True, exist_ok=True)
    try:
        morphology_profile_query.run(args)
    except ValueError as exc:
        assert "Could not resolve morphology reference bundle" in str(exc)
    else:
        raise AssertionError("Expected missing bundle resolution to fail")


def test_morphology_profile_query_neighbor_restriction_changes_output(tmp_path: Path):
    args_full = Args()
    args_full.out_dir = str(tmp_path / "morph_full")
    args_full.max_reference_neighbors = 0
    morphology_profile_query.run(args_full)

    args_limited = Args()
    args_limited.out_dir = str(tmp_path / "morph_limited")
    args_limited.max_reference_neighbors = 1
    morphology_profile_query.run(args_limited)

    full_rows = _geneset_rows(Path(args_full.out_dir) / "program=Q1__polarity=similar" / "geneset.tsv")
    limited_rows = _geneset_rows(Path(args_limited.out_dir) / "program=Q1__polarity=similar" / "geneset.tsv")
    assert [row["gene_id"] for row in full_rows] != [row["gene_id"] for row in limited_rows]


def test_morphology_profile_query_low_confidence_warning(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    args = Args()
    args.query_profiles_tsv = "tests/data/toy_morph_query_profiles_lowconf.tsv"
    args.out_dir = str(tmp_path / "morph_lowconf")
    args.query_metadata_tsv = None
    args.group_query_by = None
    args.polarity = "similar"
    result = morphology_profile_query.run(args)
    assert result["n_groups"] == 1
    captured = capsys.readouterr()
    assert "low retrieval confidence" in captured.err
    meta = json.loads((Path(args.out_dir) / "program=Q_LOW__polarity=similar" / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta["summary"]["retrieval_confidence"] == "low"


def test_morphology_profile_query_hubness_penalty_prefers_specific_neighbor(tmp_path: Path):
    query_path = tmp_path / "query.tsv"
    query_path.write_text("sample_id\tf1\tf2\nQ1\t1.0\t0.0\n", encoding="utf-8")
    ref_profiles = tmp_path / "refs.tsv"
    ref_profiles.write_text(
        "perturbation_id\tf1\tf2\n"
        "HUB\t1.0\t0.01\n"
        "SPEC\t0.8\t0.2\n",
        encoding="utf-8",
    )
    ref_meta = tmp_path / "ref_meta.tsv"
    ref_meta.write_text(
        "perturbation_id\tperturbation_type\tcompound_id\thub_score\tqc_weight\tis_control\n"
        "HUB\tcompound\tHUB\t10.0\t1.0\tfalse\n"
        "SPEC\tcompound\tSPEC\t0.1\t1.0\tfalse\n",
        encoding="utf-8",
    )
    targets = tmp_path / "targets.tsv"
    targets.write_text(
        "compound_id\tgene_symbol\tweight\n"
        "HUB\tGENE_HUB\t1.0\n"
        "SPEC\tGENE_SPEC\t1.0\n",
        encoding="utf-8",
    )

    args_none = Args()
    args_none.query_profiles_tsv = str(query_path)
    args_none.query_metadata_tsv = None
    args_none.group_query_by = None
    args_none.reference_profiles_tsv = str(ref_profiles)
    args_none.reference_metadata_tsv = str(ref_meta)
    args_none.compound_targets_tsv = str(targets)
    args_none.feature_stats_tsv = None
    args_none.feature_schema_tsv = None
    args_none.out_dir = str(tmp_path / "hub_none")
    args_none.polarity = "similar"
    args_none.top_k = 1
    args_none.gmt_topk_list = "1"
    args_none.hubness_penalty = "none"
    morphology_profile_query.run(args_none)

    args_pen = Args()
    args_pen.query_profiles_tsv = str(query_path)
    args_pen.query_metadata_tsv = None
    args_pen.group_query_by = None
    args_pen.reference_profiles_tsv = str(ref_profiles)
    args_pen.reference_metadata_tsv = str(ref_meta)
    args_pen.compound_targets_tsv = str(targets)
    args_pen.feature_stats_tsv = None
    args_pen.feature_schema_tsv = None
    args_pen.out_dir = str(tmp_path / "hub_pen")
    args_pen.polarity = "similar"
    args_pen.top_k = 1
    args_pen.gmt_topk_list = "1"
    args_pen.hubness_penalty = "inverse_rank"
    morphology_profile_query.run(args_pen)

    gene_none = _geneset_rows(Path(args_none.out_dir) / "program=Q1__polarity=similar" / "geneset.tsv")[0]["gene_id"]
    gene_pen = _geneset_rows(Path(args_pen.out_dir) / "program=Q1__polarity=similar" / "geneset.tsv")[0]["gene_id"]
    assert gene_none == "GENE_HUB"
    assert gene_pen == "GENE_SPEC"


def test_morphology_profile_query_meta_includes_specificity_fields(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "meta_specificity")
    args.polarity = "similar"
    morphology_profile_query.run(args)
    meta = json.loads((Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.meta.json").read_text(encoding="utf-8"))
    assert "neighbor_primary_target_agreement" in meta["summary"]
    assert "high_hub_mass_fraction" in meta["summary"]
    assert meta["summary"]["hubness_penalty"] == "inverse_rank"


def test_morphology_profile_query_high_retrieval_low_specificity_warns(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    query_path = tmp_path / "query.tsv"
    query_path.write_text("sample_id\tf1\tf2\nQ1\t1.0\t0.0\n", encoding="utf-8")
    ref_profiles = tmp_path / "refs.tsv"
    ref_profile_lines = ["perturbation_id\tf1\tf2"]
    ref_meta_lines = ["perturbation_id\tperturbation_type\tcompound_id\thub_score\tqc_weight\tis_control"]
    target_lines = ["compound_id\tgene_symbol\tweight"]
    for idx in range(1, 21):
        ref_profile_lines.append(f"R{idx}\t1.0\t0.0")
        ref_meta_lines.append(f"R{idx}\tcompound\tR{idx}\t0.2\t1.0\tfalse")
        target_lines.append(f"R{idx}\tGENE{idx}\t1.0")
    ref_profiles.write_text("\n".join(ref_profile_lines) + "\n", encoding="utf-8")
    ref_meta = tmp_path / "ref_meta.tsv"
    ref_meta.write_text("\n".join(ref_meta_lines) + "\n", encoding="utf-8")
    targets = tmp_path / "targets.tsv"
    targets.write_text("\n".join(target_lines) + "\n", encoding="utf-8")
    args = Args()
    args.query_profiles_tsv = str(query_path)
    args.query_metadata_tsv = None
    args.group_query_by = None
    args.reference_profiles_tsv = str(ref_profiles)
    args.reference_metadata_tsv = str(ref_meta)
    args.compound_targets_tsv = str(targets)
    args.feature_stats_tsv = None
    args.feature_schema_tsv = None
    args.out_dir = str(tmp_path / "diffuse")
    args.polarity = "similar"
    args.hubness_penalty = "none"
    args.top_k = 20
    args.gmt_topk_list = "20"
    args.gmt_min_genes = 1
    args.emit_small_gene_sets = True
    morphology_profile_query.run(args)
    captured = capsys.readouterr()
    assert "gene evidence is diffuse" in captured.err
    meta = json.loads((Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta["summary"]["retrieval_confidence"] == "high"
    assert meta["summary"]["specificity_confidence"] == "low"


def test_morphology_profile_query_opposite_experimental_warning(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    args = Args()
    args.out_dir = str(tmp_path / "opposite_warning")
    args.polarity = "opposite"
    args.min_specificity_confidence_to_emit_opposite = "low"
    morphology_profile_query.run(args)
    captured = capsys.readouterr()
    assert "opposite-polarity morphology programs are experimental" in captured.err


def test_morphology_profile_query_low_specificity_opposite_is_suppressed(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    args = Args()
    args.out_dir = str(tmp_path / "opposite_suppressed")
    args.polarity = "both"
    morphology_profile_query.run(args)
    captured = capsys.readouterr()
    assert "suppressed because specificity_confidence" in captured.err
    manifest = _manifest_rows(Path(args.out_dir) / "manifest.tsv")
    assert {row["polarity"] for row in manifest} == {"similar"}
    root_summary = json.loads((Path(args.out_dir) / "run_summary.json").read_text(encoding="utf-8"))
    assert any(row["reason"] == "low_specificity_opposite_suppressed" for row in root_summary["skipped_programs"])


def test_morphology_profile_query_can_force_opposite_emission(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "opposite_forced")
    args.polarity = "both"
    args.min_specificity_confidence_to_emit_opposite = "low"
    morphology_profile_query.run(args)
    manifest = _manifest_rows(Path(args.out_dir) / "manifest.tsv")
    assert {row["polarity"] for row in manifest} == {"similar", "opposite"}


def test_morphology_profile_query_small_gene_set_warning_includes_threshold(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    args = Args()
    args.out_dir = str(tmp_path / "small_warning")
    args.polarity = "similar"
    args.gmt_min_genes = 100
    args.gmt_max_genes = 500
    args.gmt_topk_list = "200"
    args.emit_small_gene_sets = False
    morphology_profile_query.run(args)
    captured = capsys.readouterr()
    assert "small_gene_set_skipped" in captured.err
    assert "n_genes=" in captured.err
    assert "min_required=100" in captured.err
