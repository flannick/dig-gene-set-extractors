import gzip
import json
from pathlib import Path
import tarfile
import csv

from geneset_extractors.workflows.jump_prepare_reference_bundle import run


class Args:
    profile_paths = "tests/data/toy_morph_builder_profiles.tsv"
    experimental_metadata_tsv = "tests/data/toy_morph_builder_metadata.tsv"
    compound_targets_tsv = "tests/data/toy_morph_compound_targets.tsv"
    out_dir = "tests/tmp/morph_bundle_build"
    bundle_id = "toy_jump_u2os_48h_v1"
    profile_id_column = "sample_id"
    metadata_join_id_column = "sample_id"
    perturbation_id_column = "perturbation_id"
    perturbation_type_column = "perturbation_type"
    cell_type_column = "cell_type_or_line"
    timepoint_column = "timepoint"
    gene_symbol_column = "gene_symbol"
    compound_id_column = "compound_id"
    is_control_column = "is_control"
    control_type_column = "control_type"
    cell_type_filter = "U2OS"
    timepoint_filter = "48"
    require_same_timepoint_across_modalities = True
    allow_missing_modalities = True
    allow_mixed_timepoints = False
    hubness_k = 50
    profile_kind = "normalized_feature_select_negcon_plate"
    consensus_aggregate = "median"
    profile_delimiter = "\t"
    metadata_delimiter = "\t"
    compound_targets_delimiter = "\t"
    compound_target_gene_symbol_column = "gene_symbol"
    compound_target_weight_column = "weight"
    target_annotations_tsv = None
    use_default_target_annotations = True
    target_annotations_delimiter = "\t"
    target_annotation_gene_symbol_column = "gene_symbol"


def test_jump_prepare_reference_bundle_workflow(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bundle")
    result = run(args)
    assert result["bundle_id"] == "toy_jump_u2os_48h_v1"

    out_dir = Path(args.out_dir)
    bundle_manifest = out_dir / "toy_jump_u2os_48h_v1.bundle.json"
    assert bundle_manifest.exists()
    payload = json.loads(bundle_manifest.read_text(encoding="utf-8"))
    assert payload["files"]["reference_profiles"] == "reference_profiles.tsv.gz"
    assert payload["files"]["target_annotations"] == "target_annotations.tsv.gz"
    assert payload["summary"]["n_consensus_profiles"] == 5
    assert payload["summary"]["target_annotations_source"] == "package_default"
    assert payload["summary"]["annotation_coverage"]["n_annotated_genes"] > 0
    assert "hub_score_summary" in payload["summary"]

    with gzip.open(out_dir / "reference_metadata.tsv.gz", "rt", encoding="utf-8") as fh:
        text = fh.read()
    assert "CMP_A" in text
    assert "ORF_GENE3" in text
    assert "hub_score" in text
    assert "target_family" in text
    assert (out_dir / "bundle_summary.json").exists()
    assert (out_dir / "bundle_summary.txt").exists()
    dist_dir = out_dir / "dist"
    assert (dist_dir / "dig-gene-set-extractors-toy_jump_u2os_48h_v1.tar.gz").exists()
    assert (dist_dir / "SHA256SUMS.txt").exists()
    assert (dist_dir / "distribution_manifest.json").exists()
    with tarfile.open(dist_dir / "dig-gene-set-extractors-toy_jump_u2os_48h_v1.tar.gz", "r:gz") as tf:
        members = {m.name for m in tf.getmembers()}
    assert "bundle/toy_jump_u2os_48h_v1.bundle.json" in members
    assert "bundle/reference_profiles.tsv.gz" in members
    assert "bundle/target_annotations.tsv.gz" in members
    with gzip.open(out_dir / "target_annotations.tsv.gz", "rt", encoding="utf-8") as fh:
        rows = {row["gene_symbol"]: row for row in csv.DictReader(fh, delimiter="\t")}
    assert rows["KCNN4"]["target_family"] == "Ion channel"
    assert rows["FPR1"]["target_family"] == ""


def test_jump_prepare_reference_bundle_same_timepoint_default_does_not_mix(tmp_path: Path, capsys):
    args = Args()
    args.out_dir = str(tmp_path / "bundle_same_tp")
    args.experimental_metadata_tsv = "tests/data/toy_morph_builder_metadata_mixed_timepoints.tsv"
    args.profile_paths = "tests/data/toy_morph_builder_profiles_mixed_timepoints.tsv"
    args.compound_targets_tsv = "tests/data/toy_morph_public_jump_targets.tsv"
    args.timepoint_filter = "48"
    result = run(args)
    assert result["bundle_id"] == "toy_jump_u2os_48h_v1"
    captured = capsys.readouterr()
    assert "missing modalities" in captured.err

    payload = json.loads((Path(args.out_dir) / "toy_jump_u2os_48h_v1.bundle.json").read_text(encoding="utf-8"))
    assert payload["summary"]["included_modalities"] == ["compound", "orf"]
    assert payload["summary"]["missing_modalities"] == ["crispr"]
    summary_text = (Path(args.out_dir) / "bundle_summary.txt").read_text(encoding="utf-8")
    assert "missing_modalities: ['crispr']" in summary_text
    assert "blocked_modality=crispr" in summary_text
    assert "modality=crispr included=yes" not in summary_text


def test_jump_prepare_reference_bundle_allow_mixed_timepoints(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bundle_mixed_tp")
    args.experimental_metadata_tsv = "tests/data/toy_morph_builder_metadata_mixed_timepoints.tsv"
    args.profile_paths = "tests/data/toy_morph_builder_profiles_mixed_timepoints.tsv"
    args.compound_targets_tsv = "tests/data/toy_morph_public_jump_targets.tsv"
    args.timepoint_filter = None
    args.allow_mixed_timepoints = True
    args.require_same_timepoint_across_modalities = False
    result = run(args)
    assert result["bundle_id"] == "toy_jump_u2os_48h_v1"
    payload = json.loads((Path(args.out_dir) / "toy_jump_u2os_48h_v1.bundle.json").read_text(encoding="utf-8"))
    assert payload["contexts"]["crispr"]["timepoint"] == "96"
    assert payload["contexts"]["compound"]["timepoint"] == "48"


def test_jump_prepare_reference_bundle_writes_target_annotations(tmp_path: Path):
    annotation_path = tmp_path / "target_annotations.tsv"
    annotation_path.write_text(
        "gene_symbol\ttarget_family\ttarget_class\tmechanism_label\tpathway_seed\n"
        "GENE1\tIon channel\tPotassium channel\tChannel signaling\tK_path\n"
        "GENE2\tKinase\tRTK\tNeurotrophin signaling\tN_path\n",
        encoding="utf-8",
    )
    args = Args()
    args.out_dir = str(tmp_path / "bundle_with_annotations")
    args.target_annotations_tsv = str(annotation_path)
    result = run(args)
    assert result["bundle_id"] == "toy_jump_u2os_48h_v1"
    payload = json.loads((Path(args.out_dir) / "toy_jump_u2os_48h_v1.bundle.json").read_text(encoding="utf-8"))
    assert payload["files"]["target_annotations"] == "target_annotations.tsv.gz"
    assert payload["summary"]["annotation_coverage"]["n_annotated_genes"] == 2
    assert payload["summary"]["target_annotations_source"] == "explicit"
    with gzip.open(Path(args.out_dir) / "reference_metadata.tsv.gz", "rt", encoding="utf-8") as fh:
        text = fh.read()
    assert "target_family" in text
    assert "mechanism_label" in text


def test_jump_prepare_reference_bundle_can_disable_default_target_annotations(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bundle_no_default_annotations")
    args.use_default_target_annotations = False
    result = run(args)
    assert result["bundle_id"] == "toy_jump_u2os_48h_v1"
    payload = json.loads((Path(args.out_dir) / "toy_jump_u2os_48h_v1.bundle.json").read_text(encoding="utf-8"))
    assert "target_annotations" not in payload["files"]
    assert payload["summary"]["target_annotations_source"] == "none"
    assert payload["summary"]["annotation_coverage"]["n_annotated_genes"] == 0
