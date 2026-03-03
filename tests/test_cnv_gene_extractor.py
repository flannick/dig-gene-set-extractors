import csv
import json
from pathlib import Path

import pytest

from omics2geneset.converters import cnv_gene_extractor
from omics2geneset.cli import build_parser
from omics2geneset.core.validate import validate_output_dir


def _write_toy_gtf(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                'chr1\ttest\tgene\t1000\t2000\t.\t+\t.\tgene_id "GENEA"; gene_name "GENEA"; gene_type "protein_coding";',
                'chr1\ttest\tgene\t30000000\t30001000\t.\t+\t.\tgene_id "GENEB"; gene_name "GENEB"; gene_type "protein_coding";',
                'chr1\ttest\tgene\t30002000\t30003000\t.\t+\t.\tgene_id "GENEC"; gene_name "GENEC"; gene_type "protein_coding";',
                'chr2\ttest\tgene\t5000\t6000\t.\t+\t.\tgene_id "GENED"; gene_name "GENED"; gene_type "protein_coding";',
                "",
            ]
        ),
        encoding="utf-8",
    )


class Args:
    segments_tsv = "tests/data/toy_cnv_segments.tsv"
    out_dir = "tests/tmp/cnv"
    organism = "human"
    genome_build = "hg38"
    dataset_label = "toy_cnv"
    gtf = "tests/data/toy.gtf"
    gtf_source = None
    gtf_gene_id_field = "gene_id"
    chrom_column = None
    start_column = None
    end_column = None
    amplitude_column = None
    sample_id_column = None
    segments_format = "auto"
    coord_system = "one_based_closed"
    chrom_prefix_mode = "auto"
    purity_tsv = None
    purity_sample_id_column = "sample_id"
    purity_value_column = "purity"
    use_purity_correction = False
    purity_floor = 0.1
    max_abs_amplitude = 3.0
    min_abs_amplitude = 0.10
    max_segment_length_bp = 0
    focal_length_scale_bp = 10_000_000
    focal_length_alpha = 1.0
    gene_count_penalty = "inv_sqrt"
    aggregation = "weighted_mean"
    program_preset = "default"
    program_methods = None
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
    gmt_topk_list = "3"
    gmt_mass_list = ""
    gmt_split_signed = False
    emit_small_gene_sets = True
    emit_cohort_sets = False
    cohort_score_threshold = 0.15
    cohort_min_fraction = 0.05
    cohort_min_samples = 5


def _manifest_rows(out_dir: Path) -> list[dict[str, str]]:
    with (out_dir / "manifest.tsv").open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _read_geneset_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _program_manifest_rows(out_dir: Path) -> list[dict[str, str]]:
    with (out_dir / "cnv_program_manifest.tsv").open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_cnv_gene_extractor_runs_with_focal_penalty_and_emits_gmt(tmp_path: Path):
    gtf_path = tmp_path / "toy_cnv.gtf"
    _write_toy_gtf(gtf_path)

    args = Args()
    args.gtf = str(gtf_path)
    args.out_dir = str(tmp_path / "cnv_out")
    result = cnv_gene_extractor.run(args)
    assert result["n_groups"] == 4

    out_dir = Path(args.out_dir)
    rows = _manifest_rows(out_dir)
    assert len(rows) == 4
    program_rows = _program_manifest_rows(out_dir)
    assert program_rows
    assert {
        "emitted_or_skipped",
        "n_segments_total",
        "n_segments_used",
        "n_pos_segments_used",
        "n_neg_segments_used",
        "median_seg_len_used_bp",
        "max_abs_amp_used",
        "n_genes_scored_nonzero",
        "n_genes_output",
        "primary_chrom",
    }.issubset(set(program_rows[0].keys()))
    row_lookup = {(r["sample_id"], r["program"]): r for r in rows}

    amp_s1_path = out_dir / row_lookup[("S1", "amp")]["path"] / "geneset.tsv"
    del_s1_path = out_dir / row_lookup[("S1", "del")]["path"] / "geneset.tsv"
    amp_rows = _read_geneset_rows(amp_s1_path)
    del_rows = _read_geneset_rows(del_s1_path)
    assert amp_rows
    assert del_rows
    assert amp_rows[0]["gene_id"] == "GENEA"
    assert del_rows[0]["gene_id"] == "GENED"
    assert abs(sum(float(r["weight"]) for r in amp_rows) - 1.0) < 1e-9
    assert abs(sum(float(r["weight"]) for r in del_rows) - 1.0) < 1e-9

    gmt_text = (out_dir / row_lookup[("S1", "amp")]["path"] / "genesets.gmt").read_text(encoding="utf-8")
    assert "__program=amp__topk=3" in gmt_text

    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(out_dir, schema)


def test_cnv_gene_extractor_purity_correction_increases_scores(tmp_path: Path):
    gtf_path = tmp_path / "toy_cnv.gtf"
    _write_toy_gtf(gtf_path)

    args_no = Args()
    args_no.gtf = str(gtf_path)
    args_no.program_methods = "amp"
    args_no.top_k = 3
    args_no.out_dir = str(tmp_path / "cnv_no_purity")
    cnv_gene_extractor.run(args_no)

    args_yes = Args()
    args_yes.gtf = str(gtf_path)
    args_yes.program_methods = "amp"
    args_yes.use_purity_correction = True
    args_yes.purity_tsv = "tests/data/toy_cnv_purity.tsv"
    args_yes.top_k = 3
    args_yes.out_dir = str(tmp_path / "cnv_with_purity")
    cnv_gene_extractor.run(args_yes)

    rows_no = _manifest_rows(Path(args_no.out_dir))
    rows_yes = _manifest_rows(Path(args_yes.out_dir))
    path_no = Path(args_no.out_dir) / {(r["sample_id"], r["program"]): r for r in rows_no}[("S1", "amp")]["path"] / "geneset.tsv"
    path_yes = Path(args_yes.out_dir) / {(r["sample_id"], r["program"]): r for r in rows_yes}[("S1", "amp")]["path"] / "geneset.tsv"
    score_no = next(float(r["score"]) for r in _read_geneset_rows(path_no) if r["gene_id"] == "GENEA")
    score_yes = next(float(r["score"]) for r in _read_geneset_rows(path_yes) if r["gene_id"] == "GENEA")
    assert score_yes > score_no


def test_cnv_gene_extractor_cbio_seg_auto_detection(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    gtf_path = tmp_path / "toy_cnv.gtf"
    _write_toy_gtf(gtf_path)
    segments_path = tmp_path / "cbio.seg"
    segments_path.write_text(
        "\n".join(
            [
                "ID\tchrom\tloc.start\tloc.end\tseg.mean",
                "S1\tchr1\t1100\t1300\t1.2",
                "S1\tchr2\t5100\t5300\t-1.5",
                "",
            ]
        ),
        encoding="utf-8",
    )
    args = Args()
    args.gtf = str(gtf_path)
    args.segments_tsv = str(segments_path)
    args.out_dir = str(tmp_path / "cnv_cbio")
    args.gmt_topk_list = "1"
    args.top_k = 1
    args.program_methods = "amp,del"
    args.segments_format = "auto"
    cnv_gene_extractor.run(args)
    captured = capsys.readouterr()
    assert "CNV: resolved columns:" in captured.err
    assert "sample_id=ID" in captured.err
    rows = _manifest_rows(Path(args.out_dir))
    assert ("S1", "amp") in {(r["sample_id"], r["program"]) for r in rows}
    assert ("S1", "del") in {(r["sample_id"], r["program"]) for r in rows}


def test_cnv_gene_extractor_gdc_seg_auto_detection(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    gtf_path = tmp_path / "toy_cnv.gtf"
    _write_toy_gtf(gtf_path)
    segments_path = tmp_path / "gdc.seg.tsv"
    segments_path.write_text(
        "\n".join(
            [
                "Tumor_Sample_Barcode\tChromosome\tStart\tEnd\tSegment_Mean",
                "S1\tchr1\t1100\t1300\t1.2",
                "S1\tchr2\t5100\t5300\t-1.5",
                "",
            ]
        ),
        encoding="utf-8",
    )
    args = Args()
    args.gtf = str(gtf_path)
    args.segments_tsv = str(segments_path)
    args.out_dir = str(tmp_path / "cnv_gdc")
    args.gmt_topk_list = "1"
    args.top_k = 1
    args.program_methods = "amp,del"
    args.segments_format = "auto"
    cnv_gene_extractor.run(args)
    captured = capsys.readouterr()
    assert "CNV: resolved columns:" in captured.err
    assert "sample_id=Tumor_Sample_Barcode" in captured.err
    rows = _manifest_rows(Path(args.out_dir))
    assert ("S1", "amp") in {(r["sample_id"], r["program"]) for r in rows}
    assert ("S1", "del") in {(r["sample_id"], r["program"]) for r in rows}


def test_cnv_gene_extractor_auto_ambiguous_mapping_fails(tmp_path: Path):
    gtf_path = tmp_path / "toy_cnv.gtf"
    _write_toy_gtf(gtf_path)
    segments_path = tmp_path / "ambiguous.seg.tsv"
    segments_path.write_text(
        "\n".join(
            [
                "ID\tTumor_Sample_Barcode\tchrom\tChromosome\tloc.start\tStart\tloc.end\tEnd\tseg.mean\tSegment_Mean",
                "S1\tS1\tchr1\tchr1\t1100\t1100\t1300\t1300\t1.2\t1.2",
                "",
            ]
        ),
        encoding="utf-8",
    )
    args = Args()
    args.gtf = str(gtf_path)
    args.segments_tsv = str(segments_path)
    args.out_dir = str(tmp_path / "cnv_ambiguous")
    args.segments_format = "auto"
    with pytest.raises(ValueError, match="Expected cBio SEG columns"):
        cnv_gene_extractor.run(args)


def test_cnv_gene_extractor_chr_mismatch_warning(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    gtf_path = tmp_path / "toy_cnv.gtf"
    _write_toy_gtf(gtf_path)
    segments_path = tmp_path / "mismatch.tsv"
    segments_path.write_text(
        "\n".join(
            [
                "sample_id\tchrom\tstart\tend\tsegment_mean",
                "S1\tchr1\t1100\t1300\t1.0",
                "S1\tchrUn_gl0001\t2000\t3000\t1.0",
                "S1\tchrUn_gl0002\t3000\t4000\t1.0",
                "S1\tchrUn_gl0003\t5000\t6000\t-1.0",
                "",
            ]
        ),
        encoding="utf-8",
    )

    args = Args()
    args.gtf = str(gtf_path)
    args.segments_tsv = str(segments_path)
    args.program_methods = "amp"
    args.top_k = 1
    args.gmt_topk_list = "1"
    args.out_dir = str(tmp_path / "cnv_mismatch")
    cnv_gene_extractor.run(args)
    captured = capsys.readouterr()
    assert "many segments could not be mapped to GTF chromosome names/genome build" in captured.err


def test_cnv_gene_extractor_skip_warning_includes_counts(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    gtf_path = tmp_path / "toy_cnv.gtf"
    _write_toy_gtf(gtf_path)
    segments_path = tmp_path / "only_amp.tsv"
    segments_path.write_text(
        "\n".join(
            [
                "sample_id\tChromosome\tStart\tEnd\tSegment_Mean",
                "S1\tchr1\t1100\t1300\t1.0",
                "",
            ]
        ),
        encoding="utf-8",
    )
    args = Args()
    args.gtf = str(gtf_path)
    args.segments_tsv = str(segments_path)
    args.program_methods = "amp,del"
    args.gmt_topk_list = "1"
    args.top_k = 1
    args.out_dir = str(tmp_path / "cnv_skips_warn")
    cnv_gene_extractor.run(args)
    captured = capsys.readouterr()
    assert "program=del skipped" in captured.err
    assert "counts raw=" in captured.err
    assert "after_qc=" in captured.err
    assert "overlap=" in captured.err
    program_manifest = Path(args.out_dir) / "cnv_program_manifest.tsv"
    assert program_manifest.exists()


def test_cnv_program_manifest_marks_suppressed_gmt_not_emitted(tmp_path: Path):
    gtf_path = tmp_path / "toy_cnv.gtf"
    _write_toy_gtf(gtf_path)
    args = Args()
    args.gtf = str(gtf_path)
    args.out_dir = str(tmp_path / "cnv_suppressed_gmt")
    args.program_methods = "amp"
    args.emit_small_gene_sets = False
    args.gmt_min_genes = 50
    args.gmt_max_genes = 500
    args.gmt_topk_list = "200"
    cnv_gene_extractor.run(args)

    manifest_rows = _manifest_rows(Path(args.out_dir))
    assert manifest_rows
    assert all(str(r.get("scored_ok", "")).lower() == "true" for r in manifest_rows)
    assert any(str(r.get("gmt_emitted", "")).lower() == "false" for r in manifest_rows)
    assert any(str(r.get("emit_reason", "")).strip() for r in manifest_rows)

    program_rows = _program_manifest_rows(Path(args.out_dir))
    assert any(str(r.get("status")) == "suppressed_gmt" for r in program_rows)
    assert any(str(r.get("gmt_emitted", "")).lower() == "false" for r in program_rows)
    assert all(str(r.get("status")) != "emitted_gmt" for r in program_rows if str(r.get("gmt_emitted", "")).lower() == "false")


def test_cnv_gene_extractor_writes_skipped_programs_and_git_commit(tmp_path: Path):
    gtf_path = tmp_path / "toy_cnv.gtf"
    _write_toy_gtf(gtf_path)
    segments_path = tmp_path / "only_amp.tsv"
    segments_path.write_text(
        "\n".join(
            [
                "sample_id\tChromosome\tStart\tEnd\tSegment_Mean",
                "S1\tchr1\t1100\t1300\t1.0",
                "",
            ]
        ),
        encoding="utf-8",
    )
    args = Args()
    args.gtf = str(gtf_path)
    args.segments_tsv = str(segments_path)
    args.program_methods = "amp,del"
    args.gmt_topk_list = "1"
    args.top_k = 1
    args.out_dir = str(tmp_path / "cnv_skips")
    cnv_gene_extractor.run(args)

    skipped_payload = json.loads((Path(args.out_dir) / "skipped_programs.json").read_text(encoding="utf-8"))
    assert int(skipped_payload["n_skipped"]) >= 1
    assert any(
        str(item.get("reason_code")) in {"no_positive_signal", "no_directional_input_segments"}
        for item in skipped_payload.get("items", [])
    )

    rows = _manifest_rows(Path(args.out_dir))
    amp_meta = Path(args.out_dir) / {(r["sample_id"], r["program"]): r for r in rows}[("S1", "amp")]["path"] / "geneset.meta.json"
    meta = json.loads(amp_meta.read_text(encoding="utf-8"))
    assert str(meta["converter"]["code"]["git_commit"]).strip()
    assert meta["converter"]["code"]["git_commit"] != "null"


def test_cnv_gene_extractor_tie_handling_is_deterministic(tmp_path: Path):
    gtf_path = tmp_path / "tie.gtf"
    gtf_path.write_text(
        "\n".join(
            [
                'chr1\ttest\tgene\t1000\t2000\t.\t+\t.\tgene_id "GENE_B"; gene_name "BETA"; gene_type "protein_coding";',
                'chr1\ttest\tgene\t3000\t4000\t.\t+\t.\tgene_id "GENE_A"; gene_name "ALPHA"; gene_type "protein_coding";',
                "",
            ]
        ),
        encoding="utf-8",
    )
    segments_path = tmp_path / "tie.tsv"
    segments_path.write_text(
        "\n".join(
            [
                "sample_id\tChromosome\tStart\tEnd\tSegment_Mean",
                "S1\tchr1\t1000\t4000\t1.0",
                "",
            ]
        ),
        encoding="utf-8",
    )

    def run_once(out_name: str) -> tuple[str, dict[str, object]]:
        args = Args()
        args.gtf = str(gtf_path)
        args.segments_tsv = str(segments_path)
        args.out_dir = str(tmp_path / out_name)
        args.program_methods = "amp"
        args.select = "top_k"
        args.top_k = 1
        args.gmt_topk_list = "1"
        cnv_gene_extractor.run(args)
        rows = _manifest_rows(Path(args.out_dir))
        gs = Path(args.out_dir) / rows[0]["path"] / "geneset.tsv"
        meta_path = Path(args.out_dir) / rows[0]["path"] / "geneset.meta.json"
        first_gene = _read_geneset_rows(gs)[0]["gene_id"]
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        return first_gene, meta

    g1, m1 = run_once("tie_run1")
    g2, m2 = run_once("tie_run2")
    assert g1 == g2
    tie1 = m1.get("summary", {}).get("tie_diagnostics", {})
    tie2 = m2.get("summary", {}).get("tie_diagnostics", {})
    assert tie1 == tie2
    assert bool(tie1.get("cutoff_intersects_tie_block")) is True
    assert int(tie1.get("tie_block_size", 0) or 0) >= 2


def test_cnv_gene_extractor_preflight_zero_compatibility_raises(tmp_path: Path):
    gtf_path = tmp_path / "toy_cnv.gtf"
    _write_toy_gtf(gtf_path)
    segments_path = tmp_path / "no_compat.tsv"
    segments_path.write_text(
        "\n".join(
            [
                "sample_id\tchrom\tstart\tend\tsegment_mean",
                "S1\tchrUn1\t1000\t1200\t0.8",
                "S1\tchrUn2\t2000\t2400\t-0.9",
                "",
            ]
        ),
        encoding="utf-8",
    )
    args = Args()
    args.gtf = str(gtf_path)
    args.segments_tsv = str(segments_path)
    args.program_methods = "amp"
    args.out_dir = str(tmp_path / "cnv_preflight_fail")
    with pytest.raises(ValueError, match="CNV preflight failed"):
        cnv_gene_extractor.run(args)


def test_cnv_cli_boolean_flag_ergonomics():
    parser = build_parser()
    ns_true = parser.parse_args(
        [
            "convert",
            "cnv_gene_extractor",
            "--segments_tsv",
            "tests/data/toy_cnv_segments.tsv",
            "--gtf",
            "tests/data/toy.gtf",
            "--out_dir",
            "tests/tmp/cnv_bool_true",
            "--organism",
            "human",
            "--genome_build",
            "hg38",
            "--use_purity_correction",
        ]
    )
    assert ns_true.use_purity_correction is True
    ns_false = parser.parse_args(
        [
            "convert",
            "cnv_gene_extractor",
            "--segments_tsv",
            "tests/data/toy_cnv_segments.tsv",
            "--gtf",
            "tests/data/toy.gtf",
            "--out_dir",
            "tests/tmp/cnv_bool_false",
            "--organism",
            "human",
            "--genome_build",
            "hg38",
            "--no-use-purity-correction",
        ]
    )
    assert ns_false.use_purity_correction is False
