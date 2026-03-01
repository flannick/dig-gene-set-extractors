import csv
from pathlib import Path

import pytest

from omics2geneset.converters import rna_sc_programs
from omics2geneset.core.validate import validate_output_dir


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


def test_rna_sc_programs_grouped_output_and_validation(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    args = Args()
    args.out_dir = str(tmp_path / "rna_sc_programs")
    result = rna_sc_programs.run(args)
    assert result["n_groups"] == 2

    out_dir = Path(args.out_dir)
    assert (out_dir / "manifest.tsv").exists()
    assert (out_dir / "genesets.gmt").exists()

    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    validate_output_dir(out_dir, schema)

    with (out_dir / "manifest.tsv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 2
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
