import csv
import json
from pathlib import Path

import pytest

from geneset_extractors.converters import rna_sc_programs
from geneset_extractors.core.validate import validate_output_dir
from geneset_extractors.extractors.rnaseq.sc_program_workflow import read_program_loadings
from tests.provenance_helpers import assert_manifest_has_enriched_columns


class Args:
    out_dir = "tests/tmp/rna_sc_programs"
    organism = "human"
    genome_build = "hg38"
    program_loadings_tsv = "tests/data/toy_program_loadings.tsv"
    loadings_format = "wide_genes_by_program"
    gene_id_column = "gene_id"
    program_id_column = "program_id"
    loading_column = "loading"
    transpose = False
    cnmf_gene_spectra_tsv = None
    cnmf_kind = None
    schpf_gene_scores_tsv = None
    dataset_label = "toy_sc_programs"
    signature_name = "programs"
    program_id_prefix = None
    exclude_gene_regex = None
    disable_default_excludes = False
    gtf = None
    gtf_gene_id_field = "gene_id"
    gtf_source = None
    select = "top_k"
    top_k = 200
    quantile = 0.01
    min_score = 0.0
    score_transform = "positive"
    normalize = "within_set_l1"
    emit_full = True
    emit_gmt = True
    gmt_out = None
    gmt_prefer_symbol = True
    gmt_require_symbol = False
    gmt_biotype_allowlist = "protein_coding"
    gmt_min_genes = 1
    gmt_max_genes = 10
    gmt_topk_list = "3"
    gmt_mass_list = ""
    gmt_split_signed = False
    emit_small_gene_sets = True
    provenance_overlay_json = None


def test_rna_sc_programs_grouped_output_and_validation(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    args = Args()
    args.out_dir = str(tmp_path / "rna_sc_programs")
    result = rna_sc_programs.run(args)
    assert result["n_groups"] == 2

    out_dir = Path(args.out_dir)
    assert (out_dir / "manifest.tsv").exists()
    assert (out_dir / "genesets.gmt").exists()

    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    validate_output_dir(out_dir, schema)

    with (out_dir / "manifest.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 2
    assert_manifest_has_enriched_columns(rows)
    for row in rows:
        p = out_dir / str(row["path"])
        geneset = p / "geneset.tsv"
        with geneset.open("r", encoding="utf-8") as fh:
            program_rows = list(csv.DictReader(fh, delimiter="\t"))
        assert program_rows
        weight_sum = sum(float(r["weight"]) for r in program_rows)
        assert abs(weight_sum - 1.0) < 1e-9

    gmt_lines = [line for line in (out_dir / "genesets.gmt").read_text(encoding="utf-8").splitlines() if line]
    assert gmt_lines
    for line in gmt_lines:
        name, _genes = line.split("\t")
        assert " " not in name
        assert "/" not in name

    captured = capsys.readouterr()
    assert "score_transform=positive drops negatives" in captured.err


def test_rna_sc_programs_provenance_overlay_propagates(tmp_path: Path):
    overlay_path = tmp_path / "overlay.json"
    overlay_path.write_text(
        '{"inputs":{"role:program_loadings_tsv":{"download_url":"https://example.org/programs.tsv","access_level":"public"}}}',
        encoding="utf-8",
    )
    args = Args()
    args.out_dir = str(tmp_path / "rna_sc_programs_overlay")
    args.provenance_overlay_json = str(overlay_path)
    rna_sc_programs.run(args)

    manifest_rows = list(csv.DictReader((Path(args.out_dir) / "manifest.tsv").open("r", encoding="utf-8"), delimiter="\t"))
    provenance = json.loads((Path(args.out_dir) / manifest_rows[0]["provenance_path"]).read_text(encoding="utf-8"))
    node = next(node for node in provenance["nodes"] if node.get("role") == "program_loadings_tsv")
    assert node["access"]["download_url"] == "https://example.org/programs.tsv"


def test_rna_sc_programs_signed_split_gmt(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "rna_sc_programs_signed")
    args.score_transform = "signed"
    args.gmt_split_signed = True
    args.gmt_topk_list = "2"
    result = rna_sc_programs.run(args)
    assert result["n_groups"] == 2

    out_dir = Path(args.out_dir)
    with (out_dir / "manifest.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows
    for row in rows:
        gmt_path = out_dir / str(row["path"]) / "genesets.gmt"
        text = gmt_path.read_text(encoding="utf-8")
        assert "__pos__" in text
        assert "__neg__" in text


def test_rna_sc_programs_cnmf_convenience_parser(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "rna_sc_programs_cnmf")
    args.program_loadings_tsv = None
    args.cnmf_gene_spectra_tsv = "tests/data/toy.gene_spectra_tpm.k_2.dt_0_01.txt"
    args.loadings_format = "auto"
    result = rna_sc_programs.run(args)
    assert result["n_groups"] == 2
    out_dir = Path(args.out_dir)
    assert (out_dir / "manifest.tsv").exists()
    assert (out_dir / "genesets.gmt").exists()


def test_rna_sc_programs_cnmf_programs_by_gene_orientation(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "rna_sc_programs_cnmf_programs_by_gene")
    args.program_loadings_tsv = None
    args.cnmf_gene_spectra_tsv = "tests/data/toy.gene_spectra_tpm.programs_by_gene.k_2.dt_0_01.txt"
    args.loadings_format = "auto"
    args.gmt_topk_list = "2"
    result = rna_sc_programs.run(args)
    assert result["n_groups"] == 2
    out_dir = Path(args.out_dir)
    with (out_dir / "manifest.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 2
    program_ids = {str(r["program_id"]) for r in rows}
    assert program_ids == {"1", "2"}


def test_rna_sc_programs_skips_empty_gmt_file_when_no_sets_emitted(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "rna_sc_programs_no_gmt")
    args.gmt_min_genes = 100
    args.gmt_max_genes = 500
    args.gmt_topk_list = "200"
    args.emit_small_gene_sets = False
    result = rna_sc_programs.run(args)
    assert result["n_groups"] == 2
    out_dir = Path(args.out_dir)
    with (out_dir / "manifest.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows
    for row in rows:
        gmt_path = out_dir / str(row["path"]) / "genesets.gmt"
        assert not gmt_path.exists()


def test_rna_sc_programs_schpf_convenience_parser(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "rna_sc_programs_schpf")
    args.program_loadings_tsv = None
    args.schpf_gene_scores_tsv = "tests/data/toy_program_loadings.tsv"
    args.loadings_format = "auto"
    result = rna_sc_programs.run(args)
    assert result["n_groups"] == 2
    out_dir = Path(args.out_dir)
    assert (out_dir / "manifest.tsv").exists()
    assert (out_dir / "genesets.gmt").exists()


def test_read_program_loadings_long_tidy_streaming_counts(tmp_path: Path):
    path = tmp_path / "loadings_long.tsv"
    n_rows = 0
    n_values_parsed = 0
    expected: dict[tuple[str, str], float] = {}

    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["gene_id", "program_id", "loading"])
        for i in range(3000):
            gene_id = f"G{i % 7}"
            program_id = f"P{i % 3}"
            value = float((i % 5) - 2)
            writer.writerow([gene_id, program_id, value])
            n_rows += 1
            n_values_parsed += 1
            key = (program_id, gene_id)
            expected[key] = float(expected.get(key, 0.0)) + value
        writer.writerow(["G0", "P0", "NA"])
        n_rows += 1
        writer.writerow(["", "P0", "1.0"])
        n_rows += 1
        writer.writerow(["G1", "", "2.0"])
        n_rows += 1

    programs, summary = read_program_loadings(
        path=path,
        loadings_format="long_tidy",
        gene_id_column="gene_id",
        program_id_column="program_id",
        loading_column="loading",
        transpose=False,
    )

    assert int(summary["n_rows"]) == n_rows
    assert int(summary["n_values_parsed"]) == n_values_parsed
    assert int(summary["n_values_non_numeric"]) == 1
    assert int(summary["n_programs"]) == 3
    assert int(summary["n_genes_unique"]) == 7
    assert float(programs["P0"]["G0"]) == pytest.approx(expected[("P0", "G0")])
    assert float(programs["P2"]["G6"]) == pytest.approx(expected[("P2", "G6")])


def test_rna_sc_programs_large_input_warnings(tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch):
    fake_programs = {
        "P1": {
            "GeneA": 2.0,
            "GeneB": 1.0,
            "GeneC": 0.5,
        }
    }
    fake_parse_summary = {
        "input_path": "tests/data/toy_program_loadings.tsv",
        "requested_loadings_format": "auto",
        "resolved_loadings_format": "wide_genes_by_program",
        "effective_loadings_format": "wide_genes_by_program",
        "transpose": False,
        "n_rows": 10,
        "n_programs": 250,
        "n_genes_unique": 3,
        "n_values_parsed": 6_500_000,
        "n_values_non_numeric": 0,
        "used_gene_column": "gene_id",
        "used_program_column": None,
    }

    def _fake_read_program_loadings(**_kwargs):
        return fake_programs, fake_parse_summary

    monkeypatch.setattr(rna_sc_programs, "read_program_loadings", _fake_read_program_loadings)

    args = Args()
    args.out_dir = str(tmp_path / "rna_sc_programs_large_warning")
    result = rna_sc_programs.run(args)
    assert result["n_groups"] == 1

    captured = capsys.readouterr()
    err = captured.err
    assert "large program-loading input detected" in err
    assert "very large parsed loading table detected" in err
    assert "Best practice: run factorization within cell types" in err
