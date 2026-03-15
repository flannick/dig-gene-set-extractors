import json
import hashlib
from pathlib import Path
import subprocess
import sys

from geneset_extractors.cli import build_parser


def _run(*args: str):
    return subprocess.run(
        [sys.executable, "-m", "geneset_extractors.cli", *args],
        capture_output=True,
        text=True,
        env={**__import__("os").environ, "PYTHONPATH": "src"},
    )


def test_cli_list():
    p = _run("list")
    assert p.returncode == 0
    assert "atac_bulk" in p.stdout
    assert "atac_bulk_matrix" in p.stdout
    assert "calr_ontology_mapper" in p.stdout
    assert "calr_profile_query" in p.stdout
    assert "methylation_cpg_diff" in p.stdout
    assert "methylation_dmr_regions" in p.stdout
    assert "cnv_gene_extractor" in p.stdout
    assert "drug_response_screen" in p.stdout
    assert "morphology_profile_query" in p.stdout
    assert "ptm_site_diff" in p.stdout
    assert "ptm_site_matrix" in p.stdout
    assert "rna_deg" in p.stdout
    assert "rna_deg_multi" in p.stdout
    assert "rna_sc_programs" in p.stdout
    assert "splice_event_diff" in p.stdout
    assert "splice_event_matrix" in p.stdout


def test_cli_describe():
    p = _run("describe", "atac_bulk")
    assert p.returncode == 0
    payload = json.loads(p.stdout)
    assert payload["name"] == "atac_bulk"


def test_cli_describe_bulk_matrix():
    p = _run("describe", "atac_bulk_matrix")
    assert p.returncode == 0
    payload = json.loads(p.stdout)
    assert payload["name"] == "atac_bulk_matrix"


def test_cli_describe_rna_deg():
    p = _run("describe", "rna_deg")
    assert p.returncode == 0
    payload = json.loads(p.stdout)
    assert payload["name"] == "rna_deg"
    assert any(param["name"] == "score_mode" for param in payload["parameters"])


def test_cli_describe_rna_deg_multi():
    p = _run("describe", "rna_deg_multi")
    assert p.returncode == 0
    payload = json.loads(p.stdout)
    assert payload["name"] == "rna_deg_multi"


def test_cli_describe_rna_sc_programs():
    p = _run("describe", "rna_sc_programs")
    assert p.returncode == 0
    payload = json.loads(p.stdout)
    assert payload["name"] == "rna_sc_programs"
    param_names = {str(param.get("name")) for param in payload.get("parameters", [])}
    assert "loadings_format" in param_names
    assert "score_transform" in param_names
    assert "cnmf_gene_spectra_tsv" not in param_names  # input entry, not parameter


def test_cli_describe_methylation_cpg_diff():
    p = _run("describe", "methylation_cpg_diff")
    assert p.returncode == 0
    payload = json.loads(p.stdout)
    assert payload["name"] == "methylation_cpg_diff"
    param_names = {str(param.get("name")) for param in payload.get("parameters", [])}
    assert "probe_manifest_resource_id" in param_names
    assert "delta_orientation" in param_names
    assert "exclude_gene_symbol_regex" in param_names
    assert "resources_dir" in param_names


def test_cli_describe_methylation_dmr_regions():
    p = _run("describe", "methylation_dmr_regions")
    assert p.returncode == 0
    payload = json.loads(p.stdout)
    assert payload["name"] == "methylation_dmr_regions"
    param_names = {str(param.get("name")) for param in payload.get("parameters", [])}
    assert "delta_orientation" in param_names
    assert "exclude_gene_symbol_regex" in param_names
    assert "resource_policy" in param_names


def test_cli_describe_cnv_gene_extractor():
    p = _run("describe", "cnv_gene_extractor")
    assert p.returncode == 0
    payload = json.loads(p.stdout)
    assert payload["name"] == "cnv_gene_extractor"
    param_names = {str(param.get("name")) for param in payload.get("parameters", [])}
    assert "focal_length_scale_bp" in param_names
    assert "gene_count_penalty" in param_names
    assert "use_purity_correction" in param_names


def test_cli_describe_drug_response_screen():
    p = _run("describe", "drug_response_screen")
    assert p.returncode == 0
    payload = json.loads(p.stdout)
    assert payload["name"] == "drug_response_screen"
    param_names = {str(param.get("name")) for param in payload.get("parameters", [])}
    assert "response_metric" in param_names
    assert "contrast_method" in param_names
    assert "scoring_model" in param_names
    assert "min_group_size" in param_names
    assert "target_ubiquity_penalty" in param_names
    assert "gmt_format" in param_names


def test_cli_describe_calr_ontology_mapper():
    p = _run("describe", "calr_ontology_mapper")
    assert p.returncode == 0
    payload = json.loads(p.stdout)
    assert payload["name"] == "calr_ontology_mapper"
    param_names = {str(param.get("name")) for param in payload.get("parameters", [])}
    assert "analysis_start_hour" in param_names
    assert "reference_bundle_id" in param_names
    assert "term_templates_tsv" not in param_names


def test_cli_describe_calr_profile_query():
    p = _run("describe", "calr_profile_query")
    assert p.returncode == 0
    payload = json.loads(p.stdout)
    assert payload["name"] == "calr_profile_query"
    param_names = {str(param.get("name")) for param in payload.get("parameters", [])}
    assert "similarity_metric" in param_names
    assert "reference_bundle_id" in param_names
    assert "reference_profiles_tsv" not in param_names


def test_cli_describe_morphology_profile_query():
    p = _run("describe", "morphology_profile_query")
    assert p.returncode == 0
    payload = json.loads(p.stdout)
    assert payload["name"] == "morphology_profile_query"
    param_names = {str(param.get("name")) for param in payload.get("parameters", [])}
    assert "reference_bundle_id" in param_names
    assert "similarity_metric" in param_names
    assert "polarity" in param_names


def test_cli_describe_ptm_site_diff():
    p = _run("describe", "ptm_site_diff")
    assert p.returncode == 0
    payload = json.loads(p.stdout)
    assert payload["name"] == "ptm_site_diff"
    param_names = {str(param.get("name")) for param in payload.get("parameters", [])}
    assert "protein_adjustment" in param_names
    assert "site_alias_resource_id" in param_names
    assert "gene_aggregation" in param_names
    assert "protein_adjustment_run_mode" not in param_names
    assert "emit_gene_topk_site_comparison" not in param_names


def test_cli_describe_ptm_site_matrix():
    p = _run("describe", "ptm_site_matrix")
    assert p.returncode == 0
    payload = json.loads(p.stdout)
    assert payload["name"] == "ptm_site_matrix"
    param_names = {str(param.get("name")) for param in payload.get("parameters", [])}
    assert "study_contrast" in param_names
    assert "sample_metadata_tsv" not in param_names
    assert "protein_adjustment" in param_names
    assert "protein_adjustment_run_mode" in param_names
    assert "emit_gene_topk_site_comparison" in param_names


def test_cli_describe_splice_event_diff():
    p = _run("describe", "splice_event_diff")
    assert p.returncode == 0
    payload = json.loads(p.stdout)
    assert payload["name"] == "splice_event_diff"
    param_names = {str(param.get("name")) for param in payload.get("parameters", [])}
    assert "tool_family" in param_names
    assert "impact_mode" in param_names
    assert "event_alias_resource_id" in param_names
    assert "locus_density_penalty_mode" in param_names


def test_cli_describe_splice_event_matrix():
    p = _run("describe", "splice_event_matrix")
    assert p.returncode == 0
    payload = json.loads(p.stdout)
    assert payload["name"] == "splice_event_matrix"
    param_names = {str(param.get("name")) for param in payload.get("parameters", [])}
    assert "study_contrast" in param_names
    assert "effect_metric" in param_names
    assert "sample_metadata_tsv" not in param_names
    assert "locus_density_penalty_mode" in param_names


def test_cli_workflow_ptm_prepare_public(tmp_path: Path):
    out = tmp_path / "ptm_prepare_public_cli"
    p = _run(
        "workflows",
        "ptm_prepare_public",
        "--input_mode",
        "cdap_files",
        "--ptm_report_tsv",
        "tests/data/toy_cdap_phosphosite.tmt11.tsv",
        "--protein_report_tsv",
        "tests/data/toy_cdap_proteome.tmt11.tsv",
        "--sample_design_tsv",
        "tests/data/toy_cdap.sample.txt",
        "--sample_annotations_tsv",
        "tests/data/toy_ptm_public_sample_annotations.tsv",
        "--out_dir",
        str(out),
        "--organism",
        "human",
        "--ptm_type",
        "phospho",
        "--study_id",
        "PTM_STUDY_1",
        "--study_label",
        "Toy PTM Public Study",
    )
    assert p.returncode == 0
    assert "workflow=ptm_prepare_public" in p.stderr
    assert (out / "ptm_matrix.tsv").exists()
    assert (out / "sample_metadata.tsv").exists()
    assert (out / "prepare_summary.json").exists()


def test_cli_workflow_calr_prepare_public(tmp_path: Path):
    out = tmp_path / "calr_prepare_public_cli"
    p = _run(
        "workflows",
        "calr_prepare_public",
        "--studies_tsv",
        "tests/data/toy_calr_public_studies.tsv",
        "--out_dir",
        str(out),
        "--organism",
        "mouse",
        "--bundle_id",
        "toy_calr_public_bundle_v1",
        "--write_distribution_artifact",
        "false",
    )
    assert p.returncode == 0
    assert "workflow=calr_prepare_public" in p.stderr
    assert (out / "reference_profiles.tsv").exists()
    assert (out / "reference_metadata.tsv").exists()
    assert (out / "bundle" / "toy_calr_public_bundle_v1.bundle.json").exists()


def test_cli_workflow_splice_prepare_public(tmp_path: Path):
    out = tmp_path / "splice_prepare_public_cli"
    p = _run(
        "workflows",
        "splice_prepare_public",
        "--input_mode",
        "tcga_spliceseq",
        "--psi_tsv",
        "tests/data/toy_spliceseq_public.tsv",
        "--sample_annotations_tsv",
        "tests/data/toy_spliceseq_sample_annotations.tsv",
        "--out_dir",
        str(out),
        "--organism",
        "human",
        "--genome_build",
        "hg38",
        "--study_id",
        "TCGA_TOY",
        "--study_label",
        "Toy TCGA SpliceSeq",
    )
    assert p.returncode == 0
    assert "workflow=splice_prepare_public" in p.stderr
    assert (out / "psi_matrix.tsv").exists()
    assert (out / "event_metadata.tsv").exists()
    assert (out / "prepare_summary.json").exists()


def test_cli_validate_fails_on_malformed(tmp_path: Path):
    bad = tmp_path / "bad"
    bad.mkdir()
    p = _run("validate", str(bad))
    assert p.returncode != 0


def test_cli_validate_single_good_output(tmp_path: Path):
    out = tmp_path / "bulk_cli"
    convert = _run(
        "convert",
        "atac_bulk",
        "--peaks",
        "tests/data/toy_peaks.bed",
        "--gtf",
        "tests/data/toy.gtf",
        "--out_dir",
        str(out),
        "--organism",
        "human",
        "--genome_build",
        "hg38",
        "--peak_weights_tsv",
        "tests/data/toy_peak_weights.tsv",
    )
    assert convert.returncode == 0
    assert "program_methods_active=" in convert.stderr
    assert "program_methods_skipped=" in convert.stderr
    assert "calibration_methods_active=" in convert.stderr
    assert "calibration_methods_skipped=" in convert.stderr
    validate = _run("validate", str(out))
    assert validate.returncode == 0
    assert "ok" in validate.stdout


def test_cli_validate_grouped_root(tmp_path: Path):
    out = tmp_path / "sc_grouped_cli"
    convert = _run(
        "convert",
        "atac_sc_10x",
        "--matrix_dir",
        "tests/data/toy_10x_mtx",
        "--gtf",
        "tests/data/toy.gtf",
        "--out_dir",
        str(out),
        "--organism",
        "human",
        "--genome_build",
        "hg38",
        "--groups_tsv",
        "tests/data/barcode_groups.tsv",
        "--min_cells_per_group",
        "1",
    )
    assert convert.returncode == 0
    assert "groups=" in convert.stderr
    assert "genes_per_group=" in convert.stderr
    assert "unique_genes=" in convert.stderr
    validate = _run("validate", str(out))
    assert validate.returncode == 0
    assert "n_groups=" in validate.stdout


def test_cli_convert_rna_deg(tmp_path: Path):
    out = tmp_path / "rna_deg_cli"
    convert = _run(
        "convert",
        "rna_deg",
        "--deg_tsv",
        "tests/data/toy_deg.tsv",
        "--out_dir",
        str(out),
        "--organism",
        "human",
        "--genome_build",
        "hg38",
        "--gmt_min_genes",
        "1",
        "--gmt_max_genes",
        "10",
        "--gmt_topk_list",
        "3",
        "--emit_small_gene_sets",
        "true",
    )
    assert convert.returncode == 0
    validate = _run("validate", str(out))
    assert validate.returncode == 0


def test_cli_convert_ptm_site_diff(tmp_path: Path):
    resources = tmp_path / "resources"
    resources.mkdir()
    for name in ("phosphosite_aliases_human_v1.tsv.gz", "phosphosite_ubiquity_human_v1.tsv.gz"):
        (resources / name).write_bytes((Path("tests/data") / name).read_bytes())
    out = tmp_path / "ptm_cli"
    convert = _run(
        "convert",
        "ptm_site_diff",
        "--ptm_tsv",
        "tests/data/toy_ptm_site_diff.tsv",
        "--out_dir",
        str(out),
        "--organism",
        "human",
        "--genome_build",
        "human",
        "--resources_dir",
        str(resources),
        "--gmt_min_genes",
        "1",
        "--gmt_max_genes",
        "10",
        "--gmt_topk_list",
        "2",
        "--emit_small_gene_sets",
        "true",
    )
    assert convert.returncode == 0
    validate = _run("validate", str(out))
    assert validate.returncode == 0


def test_cli_convert_ptm_site_matrix(tmp_path: Path):
    resources = tmp_path / "resources"
    resources.mkdir()
    for name in ("phosphosite_aliases_human_v1.tsv.gz", "phosphosite_ubiquity_human_v1.tsv.gz"):
        (resources / name).write_bytes((Path("tests/data") / name).read_bytes())
    out = tmp_path / "ptm_matrix_cli"
    convert = _run(
        "convert",
        "ptm_site_matrix",
        "--ptm_matrix_tsv",
        "tests/data/toy_ptm_matrix.tsv",
        "--sample_metadata_tsv",
        "tests/data/toy_ptm_sample_metadata.tsv",
        "--protein_matrix_tsv",
        "tests/data/toy_protein_matrix.tsv",
        "--study_contrast",
        "condition_within_group",
        "--out_dir",
        str(out),
        "--organism",
        "human",
        "--genome_build",
        "human",
        "--resources_dir",
        str(resources),
        "--gmt_min_genes",
        "1",
        "--gmt_max_genes",
        "20",
        "--gmt_topk_list",
        "3",
        "--emit_small_gene_sets",
        "true",
    )
    assert convert.returncode == 0
    validate = _run("validate", str(out))
    assert validate.returncode == 0


def test_cli_convert_rna_deg_multi(tmp_path: Path):
    out = tmp_path / "rna_deg_multi_cli"
    convert = _run(
        "convert",
        "rna_deg_multi",
        "--deg_tsv",
        "tests/data/toy_deg_long.tsv",
        "--comparison_column",
        "comparison_id",
        "--out_dir",
        str(out),
        "--organism",
        "human",
        "--genome_build",
        "hg38",
        "--gmt_min_genes",
        "1",
        "--gmt_max_genes",
        "10",
        "--gmt_topk_list",
        "2",
        "--emit_small_gene_sets",
        "true",
    )
    assert convert.returncode == 0
    validate = _run("validate", str(out))
    assert validate.returncode == 0


def test_cli_convert_bulk_matrix(tmp_path: Path):
    out = tmp_path / "bulk_matrix_cli"
    convert = _run(
        "convert",
        "atac_bulk_matrix",
        "--peak_matrix_tsv",
        "tests/data/toy_bulk_peak_matrix.tsv",
        "--peak_bed",
        "tests/data/toy_peaks.bed",
        "--sample_metadata_tsv",
        "tests/data/toy_bulk_sample_metadata.tsv",
        "--gtf",
        "tests/data/toy.gtf",
        "--out_dir",
        str(out),
        "--organism",
        "human",
        "--genome_build",
        "hg38",
        "--condition_a",
        "case",
        "--condition_b",
        "control",
        "--gmt_min_genes",
        "1",
        "--gmt_max_genes",
        "10",
        "--gmt_topk_list",
        "3",
        "--gmt_mass_list",
        "",
    )
    assert convert.returncode == 0
    validate = _run("validate", str(out))
    assert validate.returncode == 0


def test_cli_convert_methylation_cpg_diff(tmp_path: Path):
    out = tmp_path / "methylation_cpg_cli"
    convert = _run(
        "convert",
        "methylation_cpg_diff",
        "--cpg_tsv",
        "tests/data/toy_methylation_cpg.tsv",
        "--probe_manifest_tsv",
        "tests/data/toy_probe_manifest.tsv",
        "--gtf",
        "tests/data/toy.gtf",
        "--out_dir",
        str(out),
        "--organism",
        "human",
        "--genome_build",
        "hg38",
        "--gmt_min_genes",
        "1",
        "--gmt_max_genes",
        "10",
        "--gmt_topk_list",
        "2",
        "--emit_small_gene_sets",
        "true",
    )
    assert convert.returncode == 0
    validate = _run("validate", str(out))
    assert validate.returncode == 0


def test_cli_convert_methylation_dmr_regions(tmp_path: Path):
    out = tmp_path / "methylation_dmr_regions_cli"
    convert = _run(
        "convert",
        "methylation_dmr_regions",
        "--dmr_tsv",
        "tests/data/toy_methylation_dmr_regions.tsv",
        "--gtf",
        "tests/data/toy.gtf",
        "--out_dir",
        str(out),
        "--organism",
        "human",
        "--genome_build",
        "hg38",
        "--gmt_min_genes",
        "1",
        "--gmt_max_genes",
        "10",
        "--gmt_topk_list",
        "2",
        "--emit_small_gene_sets",
        "true",
    )
    assert convert.returncode == 0
    validate = _run("validate", str(out))
    assert validate.returncode == 0


def test_cli_convert_drug_response_screen(tmp_path: Path):
    out = tmp_path / "drug_response_cli"
    convert = _run(
        "convert",
        "drug_response_screen",
        "--response_tsv",
        "tests/data/toy_drug_response.tsv",
        "--drug_targets_tsv",
        "tests/data/toy_drug_targets.tsv",
        "--groups_tsv",
        "tests/data/toy_drug_groups.tsv",
        "--contrast_method",
        "group_vs_rest",
        "--min_group_size",
        "1",
        "--out_dir",
        str(out),
        "--organism",
        "human",
        "--genome_build",
        "hg38",
        "--gmt_min_genes",
        "1",
        "--gmt_max_genes",
        "10",
        "--gmt_topk_list",
        "2",
        "--emit_small_gene_sets",
        "true",
    )
    assert convert.returncode == 0
    validate = _run("validate", str(out))
    assert validate.returncode == 0


def test_cli_resources_list():
    p = _run("resources", "list")
    assert p.returncode == 0
    assert "encode_ccre_hg38" in p.stdout


def test_cli_resources_status(tmp_path: Path):
    p = _run("resources", "status", "--resources_dir", str(tmp_path / "res_cache"))
    assert p.returncode == 0
    assert "manual_missing" in p.stdout


def test_cli_resources_status_check_schema_known_ids():
    p = _run(
        "resources",
        "status",
        "--manifest",
        "tests/data/toy_resources_manifest.json",
        "--manifest_mode",
        "replace",
        "--resources_dir",
        "tests/data",
        "--check_schema",
    )
    assert p.returncode == 0
    rows = {}
    for line in p.stdout.splitlines():
        cols = line.split("\t")
        rows[cols[0]] = cols
    assert rows["ccre_ubiquity_hg38"][5] == "ok"
    assert rows["atac_reference_profiles_hg38"][5] == "ok"


def test_cli_resources_describe(tmp_path: Path):
    p = _run("resources", "describe", "encode_ccre_hg38", "--resources_dir", str(tmp_path / "res_cache"))
    assert p.returncode == 0
    payload = json.loads(p.stdout)
    assert payload["id"] == "encode_ccre_hg38"
    assert payload["availability"] == "manual"


def test_cli_resources_manifest_validate():
    p = _run("resources", "manifest-validate")
    assert p.returncode == 0
    assert "ok" in p.stdout


def test_cli_resources_fetch_local_manifest(tmp_path: Path):
    src_file = tmp_path / "toy_resource.tsv"
    src_file.write_text("gene\tvalue\nG1\t1\n", encoding="utf-8")
    sha = hashlib.sha256(src_file.read_bytes()).hexdigest()

    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "resources": [
                    {
                        "id": "toy_resource",
                        "stable_id": "toy-resource",
                        "version": "1",
                        "filename": "toy_resource.tsv",
                        "url": src_file.resolve().as_uri(),
                        "sha256": sha,
                        "license": "test",
                        "provider": "test",
                        "genome_build": "na",
                    }
                ],
                "presets": {"toy": ["toy_resource"]},
            }
        ),
        encoding="utf-8",
    )

    out_dir = tmp_path / "fetched"
    p = _run(
        "resources",
        "fetch",
        "--manifest",
        str(manifest),
        "--resources_dir",
        str(out_dir),
        "--preset",
        "toy",
    )
    assert p.returncode == 0
    assert "toy_resource" in p.stdout
    fetched = out_dir / "toy_resource.tsv"
    assert fetched.exists()
    assert fetched.read_text(encoding="utf-8") == src_file.read_text(encoding="utf-8")

    status = _run(
        "resources",
        "status",
        "--manifest",
        str(manifest),
        "--resources_dir",
        str(out_dir),
    )
    assert status.returncode == 0
    assert "toy_resource\tok\t" in status.stdout


def test_cli_default_gmt_bounds():
    parser = build_parser()
    args = parser.parse_args(
        [
            "convert",
            "atac_bulk",
            "--peaks",
            "tests/data/toy_peaks.bed",
            "--gtf",
            "tests/data/toy.gtf",
            "--out_dir",
            "tests/tmp/unused",
            "--organism",
            "human",
            "--genome_build",
            "hg38",
        ]
    )
    assert args.gmt_min_genes == 100
    assert args.gmt_max_genes == 500


def test_cli_resources_fetch_skips_missing_url_by_default(tmp_path: Path):
    manifest = tmp_path / "manifest_manual.json"
    manifest.write_text(
        json.dumps(
            {
                "resources": [
                    {
                        "id": "manual_resource",
                        "stable_id": "manual-resource",
                        "version": "1",
                        "filename": "manual.tsv",
                        "url": "",
                        "sha256": "",
                        "license": "test",
                        "provider": "test",
                        "genome_build": "na",
                    }
                ],
                "presets": {"manual": ["manual_resource"]},
            }
        ),
        encoding="utf-8",
    )
    out_dir = tmp_path / "fetched_manual"
    p = _run(
        "resources",
        "fetch",
        "--manifest",
        str(manifest),
        "--manifest_mode",
        "replace",
        "--resources_dir",
        str(out_dir),
        "--preset",
        "manual",
    )
    assert p.returncode == 0
    assert "manual_resource\tmanual\t" in p.stdout


def test_cli_resources_fetch_builtin_default_optional_preset_does_not_fail(tmp_path: Path):
    out_dir = tmp_path / "builtin_fetch"
    p = _run(
        "resources",
        "fetch",
        "--preset",
        "atac_default_optional",
        "--resources_dir",
        str(out_dir),
    )
    assert p.returncode == 0
    assert "\tmanual\t" in p.stdout


def test_cli_resources_status_verify_computes_checksum(tmp_path: Path):
    src_file = tmp_path / "toy_resource.tsv"
    src_file.write_text("gene\tvalue\nG1\t1\n", encoding="utf-8")
    sha = hashlib.sha256(src_file.read_bytes()).hexdigest()

    manifest = tmp_path / "manifest_verify.json"
    manifest.write_text(
        json.dumps(
            {
                "resources": [
                    {
                        "id": "toy_resource",
                        "stable_id": "toy-resource",
                        "version": "1",
                        "filename": "toy_resource.tsv",
                        "url": src_file.resolve().as_uri(),
                        "sha256": sha,
                        "license": "test",
                        "provider": "test",
                        "genome_build": "na",
                    }
                ],
                "presets": {},
            }
        ),
        encoding="utf-8",
    )
    out_dir = tmp_path / "cache"
    out_dir.mkdir()
    (out_dir / "toy_resource.tsv").write_text(src_file.read_text(encoding="utf-8"), encoding="utf-8")

    fast = _run(
        "resources",
        "status",
        "--manifest",
        str(manifest),
        "--manifest_mode",
        "replace",
        "--resources_dir",
        str(out_dir),
        "--fast",
    )
    assert fast.returncode == 0
    assert "toy_resource\tok\t" in fast.stdout

    verify = _run(
        "resources",
        "status",
        "--manifest",
        str(manifest),
        "--manifest_mode",
        "replace",
        "--resources_dir",
        str(out_dir),
        "--verify",
    )
    assert verify.returncode == 0
    assert "toy_resource\tok\t" in verify.stdout


def test_cli_resources_manifest_overlay_mode(tmp_path: Path):
    manifest = tmp_path / "overlay.json"
    manifest.write_text(
        json.dumps(
            {
                "resources": [
                    {
                        "id": "overlay_resource",
                        "stable_id": "overlay-resource",
                        "version": "1",
                        "filename": "overlay.tsv",
                        "url": "",
                        "sha256": "",
                        "license": "test",
                        "provider": "test",
                        "genome_build": "na",
                    }
                ],
                "presets": {"overlay_only": ["overlay_resource"]},
            }
        ),
        encoding="utf-8",
    )
    p = _run("resources", "list", "--manifest", str(manifest), "--manifest_mode", "overlay")
    assert p.returncode == 0
    assert "encode_ccre_hg38" in p.stdout
    assert "overlay_resource" in p.stdout


def test_cli_resources_manifest_validate_strict_fails_on_unverified_download(tmp_path: Path):
    manifest = tmp_path / "invalid_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "resources": [
                    {
                        "id": "bad_download",
                        "stable_id": "bad-download",
                        "version": "1",
                        "filename": "bad.tsv",
                        "url": "https://example.org/bad.tsv",
                        "sha256": "",
                        "license": "test",
                        "provider": "test",
                        "genome_build": "na",
                    }
                ],
                "presets": {},
            }
        ),
        encoding="utf-8",
    )
    p = _run(
        "resources",
        "manifest-validate",
        "--manifest",
        str(manifest),
        "--manifest_mode",
        "replace",
        "--strict",
    )
    assert p.returncode != 0


def test_cli_resources_status_fast_marks_unverified_without_sha(tmp_path: Path):
    src_file = tmp_path / "toy_resource.tsv"
    src_file.write_text("gene\tvalue\nG1\t1\n", encoding="utf-8")

    manifest = tmp_path / "manifest_no_sha.json"
    manifest.write_text(
        json.dumps(
            {
                "resources": [
                    {
                        "id": "toy_resource",
                        "stable_id": "toy-resource",
                        "version": "1",
                        "filename": "toy_resource.tsv",
                        "url": "",
                        "sha256": "",
                        "license": "test",
                        "provider": "test",
                        "genome_build": "na",
                    }
                ],
                "presets": {},
            }
        ),
        encoding="utf-8",
    )
    out_dir = tmp_path / "cache"
    out_dir.mkdir()
    (out_dir / "toy_resource.tsv").write_text(src_file.read_text(encoding="utf-8"), encoding="utf-8")

    fast = _run(
        "resources",
        "status",
        "--manifest",
        str(manifest),
        "--manifest_mode",
        "replace",
        "--resources_dir",
        str(out_dir),
        "--fast",
    )
    assert fast.returncode == 0
    assert "toy_resource\tok_fast\t" in fast.stdout
