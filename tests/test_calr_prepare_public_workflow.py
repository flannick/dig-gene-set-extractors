import json
from pathlib import Path
from types import SimpleNamespace

from geneset_extractors.converters import calr_profile_query
from geneset_extractors.core.validate import validate_output_dir
from geneset_extractors.workflows.calr_prepare_public import run as run_calr_prepare_public


class QueryArgs:
    calr_data_csv = "tests/data/toy_calr_data.csv"
    session_csv = "tests/data/toy_calr_session.csv"
    exclusions_tsv = None
    out_dir = ""
    organism = "mouse"
    genome_build = "mm39"
    dataset_label = "toy_query"
    signature_name = "toy_query_sig"
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
    reference_profiles_tsv = None
    reference_metadata_tsv = None
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


def test_calr_prepare_public_builds_reference_tables_and_bundle(tmp_path: Path):
    out_dir = tmp_path / "calr_public"
    result = run_calr_prepare_public(
        SimpleNamespace(
            studies_tsv="tests/data/toy_calr_public_studies.tsv",
            out_dir=str(out_dir),
            organism="mouse",
            bundle_id="toy_calr_public_bundle_v1",
            build_bundle=True,
            write_distribution_artifact=False,
            term_templates_tsv=None,
            phenotype_gene_edges_tsv=None,
            term_hierarchy_tsv=None,
            include_packaged_term_hierarchy=True,
            exploratory_without_session=True,
            min_group_size=2,
            mass_covariate=None,
            analysis_start_hour=None,
            analysis_end_hour=None,
            photoperiod_lights_on_hour=None,
            photoperiod_hours_light=12.0,
        )
    )
    assert result["n_profiles"] == 2
    assert (out_dir / "reference_profiles.tsv").exists()
    assert (out_dir / "reference_metadata.tsv").exists()
    assert (out_dir / "feature_schema.tsv").exists()
    assert (out_dir / "feature_stats.tsv").exists()
    assert (out_dir / "resolved_studies.tsv").exists()
    assert (out_dir / "prepare_summary.json").exists()
    assert (out_dir / "bundle" / "toy_calr_public_bundle_v1.bundle.json").exists()

    summary = json.loads((out_dir / "prepare_summary.json").read_text(encoding="utf-8"))
    assert summary["n_reference_profiles"] == 2
    assert summary["n_explicit_session_studies"] == 2

    metadata_lines = (out_dir / "reference_metadata.tsv").read_text(encoding="utf-8").splitlines()
    assert any("wide_membership" in line for line in metadata_lines)
    assert any("Ppargc1a" in line for line in metadata_lines)


def test_calr_prepare_public_bundle_runs_profile_query(tmp_path: Path):
    public_out = tmp_path / "calr_public"
    run_calr_prepare_public(
        SimpleNamespace(
            studies_tsv="tests/data/toy_calr_public_studies.tsv",
            out_dir=str(public_out),
            organism="mouse",
            bundle_id="toy_calr_public_bundle_v1",
            build_bundle=True,
            write_distribution_artifact=False,
            term_templates_tsv=None,
            phenotype_gene_edges_tsv=None,
            term_hierarchy_tsv=None,
            include_packaged_term_hierarchy=True,
            exploratory_without_session=True,
            min_group_size=2,
            mass_covariate=None,
            analysis_start_hour=None,
            analysis_end_hour=None,
            photoperiod_lights_on_hour=None,
            photoperiod_hours_light=12.0,
        )
    )

    args = QueryArgs()
    args.out_dir = str(tmp_path / "calr_profile_from_public")
    args.reference_bundle_id = "toy_calr_public_bundle_v1"
    args.resources_dir = str(public_out / "bundle")
    result = calr_profile_query.run(args)
    assert result["n_groups"] > 0
    validate_output_dir(Path(args.out_dir), Path("src/geneset_extractors/schemas/geneset_metadata.schema.json"))
