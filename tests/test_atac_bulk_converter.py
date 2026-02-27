import csv
import json
from pathlib import Path

import pytest

from omics2geneset.converters import atac_bulk
from omics2geneset.core.validate import validate_output_dir


class Args:
    peaks = "tests/data/toy_peaks.bed"
    gtf = "tests/data/toy.gtf"
    out_dir = "tests/tmp/bulk"
    organism = "human"
    genome_build = "hg38"
    peaks_weight_column = 5
    peak_weights_tsv = "tests/data/toy_peak_weights.tsv"
    link_method = "all"
    promoter_upstream_bp = 2000
    promoter_downstream_bp = 500
    max_distance_bp = None
    decay_length_bp = 50000
    max_genes_per_peak = 5
    peak_weight_transform = "abs"
    normalize = "within_set_l1"
    program_preset = "connectable"
    program_methods = None
    calibration_methods = "auto_prefer_ref_ubiquity_else_none"
    resources_manifest = None
    resources_dir = None
    use_reference_bundle = True
    resource_policy = "skip"
    ref_ubiquity_resource_id = None
    atlas_resource_id = None
    atlas_metric = "zscore"
    atlas_eps = 1e-6
    atlas_min_raw_quantile = 0.95
    atlas_use_log1p = True
    select = "top_k"
    top_k = 200
    quantile = 0.01
    min_score = 0.0
    emit_full = True
    emit_gmt = True
    gmt_out = None
    gmt_prefer_symbol = True
    gmt_require_symbol = True
    gmt_biotype_allowlist = "protein_coding"
    gmt_min_genes = 100
    gmt_max_genes = 500
    gmt_topk_list = "200"
    gmt_mass_list = ""
    gmt_split_signed = False
    emit_small_gene_sets = False
    qc_marker_genes_tsv = None
    region_gene_links_tsv = None
    gtf_source = "toy"
    dataset_label = None


def test_bulk_converter_end_to_end(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bulk")
    args.top_k = 1
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    atac_bulk.run(args)

    geneset = Path(args.out_dir) / "geneset.tsv"
    geneset_full = Path(args.out_dir) / "geneset.full.tsv"
    meta = Path(args.out_dir) / "geneset.meta.json"
    run_summary_json = Path(args.out_dir) / "run_summary.json"
    run_summary_txt = Path(args.out_dir) / "run_summary.txt"
    assert geneset.exists()
    assert geneset_full.exists()
    assert meta.exists()
    assert run_summary_json.exists()
    assert run_summary_txt.exists()

    with geneset.open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1
    assert abs(sum(float(r["weight"]) for r in rows) - 1.0) < 1e-9

    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(Path(args.out_dir), schema)


def test_bulk_optional_marker_qc_written(tmp_path: Path):
    markers_path = tmp_path / "markers.tsv"
    markers_path.write_text("gene_symbol\nGENE1\nGENE2\n", encoding="utf-8")

    args = Args()
    args.out_dir = str(tmp_path / "bulk_marker_qc")
    args.top_k = 2
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    args.qc_marker_genes_tsv = str(markers_path)
    atac_bulk.run(args)

    meta = json.loads((Path(args.out_dir) / "geneset.meta.json").read_text(encoding="utf-8"))
    assert "marker_qc" in meta["summary"]
    assert meta["summary"]["marker_qc"]["n_markers_provided"] == 2

    run_summary = json.loads((Path(args.out_dir) / "run_summary.json").read_text(encoding="utf-8"))
    assert "marker_qc" in run_summary


def test_bulk_default_connectable_emits_exactly_two_sets(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bulk_connectable")
    args.resources_manifest = "tests/data/toy_resources_manifest.json"
    args.resources_dir = "tests/data"
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    atac_bulk.run(args)

    lines = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == 2
    text = "\n".join(lines)
    assert "__program=linked_activity__calibration_method=ref_ubiquity_penalty__link_method=nearest_tss__topk=3" in text
    assert "__program=distal_activity__calibration_method=ref_ubiquity_penalty__link_method=distance_decay__topk=3" in text
    assert "__program=promoter_activity__" not in text
    assert "__program=enhancer_bias__" not in text
    assert "__calibration_method=atlas_residual__" not in text
    assert "__link_method=promoter_overlap__" not in text


def test_bulk_default_auto_falls_back_to_none_without_resources(tmp_path: Path, capsys):
    args = Args()
    args.out_dir = str(tmp_path / "bulk_connectable_no_resources")
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    atac_bulk.run(args)
    captured = capsys.readouterr()
    assert "auto_prefer_ref_ubiquity_else_none" in captured.err
    assert "falling back to calibration_method=none" in captured.err

    text = (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8")
    assert "__calibration_method=none__" in text
    assert "__program=linked_activity__calibration_method=none__link_method=nearest_tss__topk=3" in text
    assert "__program=distal_activity__calibration_method=none__link_method=distance_decay__topk=3" in text


def test_bulk_use_reference_bundle_false_does_not_crash(tmp_path: Path, capsys):
    args = Args()
    args.out_dir = str(tmp_path / "bulk_no_bundle")
    args.use_reference_bundle = False
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    atac_bulk.run(args)
    captured = capsys.readouterr()
    assert "auto_prefer_ref_ubiquity_else_none" in captured.err
    meta = json.loads((Path(args.out_dir) / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta["converter"]["parameters"]["use_reference_bundle"] is False


def test_bulk_ref_ubiquity_changes_distal_ranking(tmp_path: Path):
    peaks_path = tmp_path / "peaks.bed"
    peaks_path.write_text(
        "\n".join(
            [
                "chr1\t2500\t2600\tp1\t5.0",
                "chr1\t9000\t9100\tp2\t1.0",
                "",
            ]
        ),
        encoding="utf-8",
    )
    gtf_path = tmp_path / "genes.gtf"
    gtf_path.write_text(
        "\n".join(
            [
                'chr1\ttest\tgene\t1000\t1500\t.\t+\t.\tgene_id "GENE_A"; gene_name "GA";',
                'chr1\ttest\tgene\t12000\t12500\t.\t+\t.\tgene_id "GENE_B"; gene_name "GB";',
                "",
            ]
        ),
        encoding="utf-8",
    )
    ref_path = tmp_path / "ref.tsv"
    ref_path.write_text(
        "\n".join(
            [
                "chrom\tstart\tend\tidf_ref",
                "chr1\t2400\t2700\t0.01",
                "chr1\t8900\t9200\t100.0",
                "",
            ]
        ),
        encoding="utf-8",
    )
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "resources": [
                    {
                        "id": "toy_ref",
                        "description": "toy",
                        "filename": ref_path.name,
                        "url": "",
                        "sha256": "",
                        "provider": "test",
                        "stable_id": "toy-ref",
                        "version": "1",
                        "license": "test",
                        "genome_build": "hg38",
                    }
                ],
                "presets": {},
            }
        ),
        encoding="utf-8",
    )

    args_none = Args()
    args_none.peaks = str(peaks_path)
    args_none.gtf = str(gtf_path)
    args_none.peak_weights_tsv = None
    args_none.out_dir = str(tmp_path / "none")
    args_none.program_preset = "none"
    args_none.program_methods = "distal_activity"
    args_none.link_method = "distance_decay"
    args_none.calibration_methods = "none"
    args_none.gmt_min_genes = 1
    args_none.gmt_max_genes = 10
    args_none.gmt_topk_list = "1"
    atac_bulk.run(args_none)

    args_ref = Args()
    args_ref.peaks = str(peaks_path)
    args_ref.gtf = str(gtf_path)
    args_ref.peak_weights_tsv = None
    args_ref.out_dir = str(tmp_path / "ref")
    args_ref.program_preset = "none"
    args_ref.program_methods = "distal_activity"
    args_ref.link_method = "distance_decay"
    args_ref.calibration_methods = "ref_ubiquity_penalty"
    args_ref.resources_manifest = str(manifest_path)
    args_ref.resources_dir = str(tmp_path)
    args_ref.ref_ubiquity_resource_id = "toy_ref"
    args_ref.gmt_min_genes = 1
    args_ref.gmt_max_genes = 10
    args_ref.gmt_topk_list = "1"
    atac_bulk.run(args_ref)

    none_line = (Path(args_none.out_dir) / "genesets.gmt").read_text(encoding="utf-8").strip().splitlines()[0]
    ref_line = (Path(args_ref.out_dir) / "genesets.gmt").read_text(encoding="utf-8").strip().splitlines()[0]
    none_gene = none_line.split("\t", 1)[1]
    ref_gene = ref_line.split("\t", 1)[1]
    assert none_gene != ref_gene


def test_bulk_ref_ubiquity_applies_idf_once(tmp_path: Path):
    peaks_path = tmp_path / "peaks.bed"
    peaks_path.write_text(
        "\n".join(
            [
                "chr1\t1000\t1100\tp1\t2.0",
                "chr1\t2000\t2100\tp2\t3.0",
                "",
            ]
        ),
        encoding="utf-8",
    )
    gtf_path = tmp_path / "genes.gtf"
    gtf_path.write_text(
        "\n".join(
            [
                'chr1\ttest\tgene\t1450\t1550\t.\t+\t.\tgene_id "GENE1"; gene_name "G1";',
                "",
            ]
        ),
        encoding="utf-8",
    )
    ref_path = tmp_path / "ref.tsv"
    ref_path.write_text(
        "\n".join(
            [
                "chrom\tstart\tend\tidf_ref",
                "chr1\t900\t1200\t10.0",
                "chr1\t1900\t2200\t1.0",
                "",
            ]
        ),
        encoding="utf-8",
    )
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "resources": [
                    {
                        "id": "toy_ref",
                        "description": "toy",
                        "filename": ref_path.name,
                        "url": "",
                        "sha256": "",
                        "provider": "test",
                        "stable_id": "toy-ref",
                        "version": "1",
                        "license": "test",
                        "genome_build": "hg38",
                    }
                ],
                "presets": {},
            }
        ),
        encoding="utf-8",
    )

    args = Args()
    args.peaks = str(peaks_path)
    args.gtf = str(gtf_path)
    args.peak_weights_tsv = None
    args.out_dir = str(tmp_path / "ref_once")
    args.program_preset = "none"
    args.program_methods = "linked_activity"
    args.link_method = "nearest_tss"
    args.calibration_methods = "ref_ubiquity_penalty"
    args.resources_manifest = str(manifest_path)
    args.resources_dir = str(tmp_path)
    args.ref_ubiquity_resource_id = "toy_ref"
    args.top_k = 1
    args.emit_gmt = False
    atac_bulk.run(args)

    with (Path(args.out_dir) / "geneset.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["gene_id"] == "GENE1"
    assert float(rows[0]["score"]) == pytest.approx(23.0, rel=0.0, abs=1e-9)


def test_bulk_atlas_missing_score_definition_warns_and_skips(tmp_path: Path, capsys):
    atlas_path = tmp_path / "atlas.tsv"
    atlas_path.write_text(
        "\n".join(
            [
                "score_definition\tgene_id\tmedian_score\tmad_score",
                "program=linked_activity__link_method=nearest_tss\tGENE1\t1.0\t0.1",
                "program=linked_activity__link_method=nearest_tss\tGENE2\t1.0\t0.1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "resources": [
                    {
                        "id": "atlas_custom",
                        "description": "toy atlas",
                        "filename": atlas_path.name,
                        "url": "",
                        "sha256": "",
                        "provider": "test",
                        "stable_id": "toy-atlas",
                        "version": "1",
                        "license": "test",
                        "genome_build": "hg38",
                    }
                ],
                "presets": {},
            }
        ),
        encoding="utf-8",
    )

    args = Args()
    args.out_dir = str(tmp_path / "atlas_missing")
    args.calibration_methods = "atlas_residual"
    args.resources_manifest = str(manifest_path)
    args.resources_dir = str(tmp_path)
    args.atlas_resource_id = "atlas_custom"
    args.gmt_min_genes = 1
    args.gmt_max_genes = 10
    args.gmt_topk_list = "3"
    atac_bulk.run(args)
    captured = capsys.readouterr()
    assert "atlas_residual baseline missing for score definitions" in captured.err


def test_bulk_external_requires_region_gene_links(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "bulk_external_missing")
    args.link_method = "external"
    args.region_gene_links_tsv = None
    with pytest.raises(ValueError, match="region_gene_links_tsv"):
        atac_bulk.run(args)


def test_bulk_skips_small_gmt_sets_by_default(tmp_path: Path, capsys):
    args = Args()
    args.out_dir = str(tmp_path / "bulk_small_skipped")
    args.peak_weights_tsv = None
    args.gmt_min_genes = 5
    args.gmt_max_genes = 10
    args.gmt_topk_list = "5"
    args.calibration_methods = "none"
    atac_bulk.run(args)
    captured = capsys.readouterr()
    assert "skipped GMT output" in captured.err
    gmt_path = Path(args.out_dir) / "genesets.gmt"
    assert gmt_path.exists()
    assert gmt_path.read_text(encoding="utf-8").strip() == ""
    meta = json.loads((Path(args.out_dir) / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta["gmt"]["skipped_outputs"]


def test_bulk_emit_small_gene_sets_override(tmp_path: Path, capsys):
    args = Args()
    args.out_dir = str(tmp_path / "bulk_small_emitted")
    args.peak_weights_tsv = None
    args.gmt_min_genes = 5
    args.gmt_max_genes = 10
    args.gmt_topk_list = "5"
    args.emit_small_gene_sets = True
    args.calibration_methods = "none"
    atac_bulk.run(args)
    captured = capsys.readouterr()
    assert "emitted small GMT output" in captured.err
    lines = [x for x in (Path(args.out_dir) / "genesets.gmt").read_text(encoding="utf-8").splitlines() if x.strip()]
    assert lines
