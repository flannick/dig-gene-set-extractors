import json
from pathlib import Path
import subprocess
import sys


def _run(*args: str):
    return subprocess.run(
        [sys.executable, "-m", "omics2geneset.cli", *args],
        capture_output=True,
        text=True,
        env={**__import__("os").environ, "PYTHONPATH": "src"},
    )


def test_cli_list():
    p = _run("list")
    assert p.returncode == 0
    assert "atac_bulk" in p.stdout


def test_cli_describe():
    p = _run("describe", "atac_bulk")
    assert p.returncode == 0
    payload = json.loads(p.stdout)
    assert payload["name"] == "atac_bulk"


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
    )
    assert convert.returncode == 0
    assert "groups=" in convert.stderr
    assert "genes_per_group=" in convert.stderr
    assert "unique_genes=" in convert.stderr
    validate = _run("validate", str(out))
    assert validate.returncode == 0
    assert "n_groups=" in validate.stdout
