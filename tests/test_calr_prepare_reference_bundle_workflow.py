import json
from pathlib import Path

from geneset_extractors.workflows.calr_prepare_reference_bundle import run


def test_calr_prepare_reference_bundle_workflow(tmp_path: Path):
    args = type("Args", (), {
        "reference_profiles_tsv": "tests/data/calorimetry_bundle/toy_reference_profiles.tsv",
        "reference_metadata_tsv": "tests/data/calorimetry_bundle/toy_reference_metadata.tsv",
        "feature_schema_tsv": None,
        "feature_stats_tsv": None,
        "term_templates_tsv": None,
        "phenotype_gene_edges_tsv": None,
        "term_hierarchy_tsv": None,
        "include_packaged_term_hierarchy": True,
        "out_dir": str(tmp_path / "calr_bundle"),
        "organism": "mouse",
        "bundle_id": "toy_calr_bundle_v1",
        "write_distribution_artifact": True,
        "distribution_dir": None,
    })()
    result = run(args)
    out_dir = Path(args.out_dir)
    assert result["bundle_id"] == "toy_calr_bundle_v1"
    assert (out_dir / "reference_profiles.tsv.gz").exists()
    assert (out_dir / "reference_metadata.tsv.gz").exists()
    assert (out_dir / "feature_schema.tsv.gz").exists()
    assert (out_dir / "feature_stats.tsv.gz").exists()
    assert (out_dir / "term_templates.tsv.gz").exists()
    assert (out_dir / "phenotype_gene_edges.tsv.gz").exists()
    assert (out_dir / "term_hierarchy.tsv.gz").exists()
    assert (out_dir / "toy_calr_bundle_v1.bundle.json").exists()
    assert (out_dir / "bundle_summary.json").exists()
    assert (out_dir / "bundle_summary.txt").exists()
    assert (out_dir / "dist" / "SHA256SUMS.txt").exists()
    assert result["tarball"].endswith(".tar.gz")

    summary = json.loads((out_dir / "bundle_summary.json").read_text(encoding="utf-8"))
    assert summary["term_templates_source"] == "packaged_default"
    assert summary["n_reference_profiles"] == 4
