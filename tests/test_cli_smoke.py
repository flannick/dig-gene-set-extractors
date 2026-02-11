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
