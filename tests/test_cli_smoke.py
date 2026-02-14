import json
import hashlib
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


def test_cli_resources_list():
    p = _run("resources", "list")
    assert p.returncode == 0
    assert "encode_ccre_hg38" in p.stdout


def test_cli_resources_status(tmp_path: Path):
    p = _run("resources", "status", "--resources_dir", str(tmp_path / "res_cache"))
    assert p.returncode == 0
    assert "missing" in p.stdout


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
                        "filename": "toy_resource.tsv",
                        "url": src_file.resolve().as_uri(),
                        "sha256": sha,
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
