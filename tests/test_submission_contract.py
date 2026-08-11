from __future__ import annotations

import json
import os
import subprocess
import sys
import tomllib
from pathlib import Path

from geneset_extractors import submission
from geneset_extractors.cli import main


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, "-m", "geneset_extractors.cli", *args], capture_output=True, text=True, env={**os.environ, "PYTHONPATH": "src"})


def test_registered_converter_resolution_and_smoke() -> None:
    result = submission.validate_contract("rna_deg")
    assert result["ok"]
    assert result["checks"]["registered"]
    assert result["checks"]["module_importable"]
    assert result["checks"]["smoke_test"]["ok"]


def test_registered_workflow_resolution() -> None:
    result = submission.validate_contract("rna_de_prepare")
    assert result["ok"]
    assert result["contract"]["kind"] == "workflow"
    assert result["checks"]["registered"]


def test_unknown_identifier_reports_failure() -> None:
    result = submission.validate_contract("not_a_real_submission_id")
    assert not result["ok"]
    assert result["checks"]["registered"] is False


def test_import_failure_is_reported(monkeypatch) -> None:
    original = submission.describe_contract
    contract = original("rna_deg")
    contract["source_module"] = "geneset_extractors.not_importable"
    monkeypatch.setattr(submission, "describe_contract", lambda _: contract)
    result = submission.validate_contract("rna_deg")
    assert not result["ok"]
    assert result["checks"]["module_importable"] is False
    assert "ModuleNotFoundError" in result["checks"]["import_error"]


def test_submission_cli_list_describe_and_validate() -> None:
    listed = _run("submission", "list")
    assert listed.returncode == 0
    assert any(row["identifier"] == "rna_deg" for row in json.loads(listed.stdout))
    described = _run("submission", "describe", "rna_deg")
    assert described.returncode == 0
    assert json.loads(described.stdout)["cli_command"] == ["convert", "rna_deg"]
    validated = _run("submission", "validate", "rna_deg")
    assert validated.returncode == 0
    assert json.loads(validated.stdout)["ok"] is True


def test_existing_cli_and_aliases_remain_available() -> None:
    assert main(["describe", "rna_deg"]) == 0
    with (Path(__file__).resolve().parents[1] / "pyproject.toml").open("rb") as handle:
        scripts = tomllib.load(handle)["project"]["scripts"]
    assert scripts["geneset-extractors"] == "geneset_extractors.cli:main"
    assert scripts["geneset_extractors"] == "geneset_extractors.cli:main"
    assert scripts["omics2geneset"] == "geneset_extractors.cli:main"
