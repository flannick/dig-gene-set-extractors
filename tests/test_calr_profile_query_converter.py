import csv
import json
from pathlib import Path

from geneset_extractors.converters import calr_profile_query
from geneset_extractors.core.validate import validate_output_dir
from geneset_extractors.workflows.calr_prepare_reference_bundle import run as run_calr_prepare_reference_bundle
from tests.provenance_helpers import (
    assert_node_has_structured_resource_metadata,
    file_node_for_role,
    load_provenance,
)


class Args:
    calr_data_csv = "tests/data/toy_calr_data.csv"
    session_csv = "tests/data/toy_calr_session.csv"
    exclusions_tsv = "tests/data/toy_calr_exclusions.tsv"
    out_dir = "tests/tmp/calr_profile"
    organism = "mouse"
    genome_build = "mm39"
    dataset_label = "toy_calr"
    signature_name = "toy_calr_profile"
    analysis_start_hour = None
    analysis_end_hour = None
    photoperiod_lights_on_hour = None
    photoperiod_hours_light = 12.0
    exploratory_without_session = True
    mass_covariate = None
    min_group_size = 2
    resources_manifest = None
    resources_dir = None
    resource_policy = "skip"
    reference_bundle_id = None
    reference_profiles_tsv = "tests/data/calorimetry_bundle/toy_reference_profiles.tsv"
    reference_metadata_tsv = "tests/data/calorimetry_bundle/toy_reference_metadata.tsv"
    feature_schema_tsv = None
    feature_stats_tsv = None
    similarity_metric = "cosine"
    similarity_floor = 0.0
    similarity_power = 1.0
    hubness_penalty = "inverse_linear"
    provenance_mismatch_penalty = 0.15
    select = "top_k"
    top_k = 10
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
    gmt_max_genes = 20
    gmt_topk_list = "5"
    gmt_mass_list = ""
    gmt_split_signed = False
    gmt_format = "classic"
    emit_small_gene_sets = True
    output_gene_species = "human"
    ortholog_policy = "unique_only"
    mouse_human_orthologs_tsv = None
    provenance_overlay_json = None


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_calr_profile_query_explicit_reference_files(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "calr_profile")
    result = calr_profile_query.run(args)
    assert result["n_groups"] > 0

    out_dir = Path(args.out_dir)
    manifest_rows = _read_rows(out_dir / "manifest.tsv")
    assert manifest_rows
    core_path = out_dir / "program=thermogenesis__mode=core__contrast=KO" / "geneset.tsv"
    rows = _read_rows(core_path)
    assert rows
    assert rows[0]["gene_symbol"] in {"UCP1", "MLXIPL", "CLOCK", "NPY"}
    meta = json.loads((out_dir / "program=thermogenesis__mode=core__contrast=KO" / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta["summary"]["retrieval_confidence"] in {"medium", "high", "low"}
    assert meta["summary"]["output_gene_species"] == "human"
    validate_output_dir(out_dir, Path("src/geneset_extractors/schemas/geneset_metadata.schema.json"))


def test_calr_profile_query_bundle_mode(tmp_path: Path):
    bundle_args = type("BundleArgs", (), {
        "reference_profiles_tsv": "tests/data/calorimetry_bundle/toy_reference_profiles.tsv",
        "reference_metadata_tsv": "tests/data/calorimetry_bundle/toy_reference_metadata.tsv",
        "feature_schema_tsv": None,
        "feature_stats_tsv": None,
        "term_templates_tsv": None,
        "phenotype_gene_edges_tsv": None,
        "term_hierarchy_tsv": None,
        "include_packaged_term_hierarchy": True,
        "out_dir": str(tmp_path / "bundle"),
        "organism": "mouse",
        "bundle_id": "toy_calr_bundle_v1",
        "write_distribution_artifact": True,
        "distribution_dir": None,
    })()
    run_calr_prepare_reference_bundle(bundle_args)

    args = Args()
    args.out_dir = str(tmp_path / "profile_bundle")
    args.reference_profiles_tsv = None
    args.reference_metadata_tsv = None
    args.resources_dir = str(bundle_args.out_dir)
    args.reference_bundle_id = "toy_calr_bundle_v1"
    result = calr_profile_query.run(args)
    assert result["n_groups"] > 0
    meta = json.loads((Path(args.out_dir) / "program=global__mode=core__contrast=KO" / "geneset.meta.json").read_text(encoding="utf-8"))
    resources = meta["summary"]["resources"]
    assert resources is not None
    assert any(row["id"] == "toy_calr_bundle_v1" for row in resources["used"])
    assert any(row["id"] == "calorimetry_mouse_human_orthologs_v1" for row in resources["used"])
    assert meta["summary"]["output_gene_species"] == "human"
    provenance = load_provenance(Path(args.out_dir) / "program=global__mode=core__contrast=KO")
    node = file_node_for_role(provenance, "reference_profiles_tsv")
    assert_node_has_structured_resource_metadata(node)


def test_calr_profile_query_low_confidence_warning(tmp_path: Path, capsys):
    ref_profiles = tmp_path / "refs.tsv"
    ref_profiles.write_text(
        "perturbation_id\tvo2_mean_selected\tee_mean_selected\trer_mean_selected\tfeed_auc_selected\txytot_mean_selected\tbody.temp_mean_selected\txytot_amp1\n"
        "GENE_A\t0.05\t0.03\t0.02\t0.01\t0.02\t0.01\t0.02\n"
        "GENE_B\t0.02\t0.01\t0.03\t0.02\t0.01\t0.02\t0.01\n",
        encoding="utf-8",
    )
    ref_meta = tmp_path / "ref_meta.tsv"
    ref_meta.write_text(
        "perturbation_id\tgene_id\tgene_symbol\tqc_weight\thub_score\tmass_covariate\tacclimation_state\tambient_temperature\n"
        "GENE_A\tGENE_A\tGENE_A\t1.0\t0.1\tlean_mass\tpost_window_selected\t22\n"
        "GENE_B\tGENE_B\tGENE_B\t1.0\t0.1\tlean_mass\tpost_window_selected\t22\n",
        encoding="utf-8",
    )
    args = Args()
    args.out_dir = str(tmp_path / "calr_profile_low")
    args.reference_profiles_tsv = str(ref_profiles)
    args.reference_metadata_tsv = str(ref_meta)
    args.output_gene_species = "source"
    result = calr_profile_query.run(args)
    assert result["n_groups"] > 0
    captured = capsys.readouterr()
    assert "low retrieval confidence" in captured.err
    manifest_rows = _read_rows(Path(args.out_dir) / "manifest.tsv")
    low_row = manifest_rows[0]
    summary = json.loads((Path(args.out_dir) / low_row["path"] / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["retrieval_confidence"] == "low"
