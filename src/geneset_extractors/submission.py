"""Reusable DIG-side contract inspection for wrapper submission validation.

This module deliberately exposes existing CLI/registry behavior; it does not
contain wrapper orchestration, publishing, or a second workflow registry.
"""
from __future__ import annotations

import argparse
import importlib
import tempfile
from pathlib import Path
from typing import Any

from geneset_extractors.core.validate import validate_output_dir
from geneset_extractors.registry import CONVERTERS, get_converter, get_converter_spec


def _subparser_choices(parser: argparse.ArgumentParser, destination: str) -> dict[str, argparse.ArgumentParser]:
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction) and action.dest == destination:
            return dict(action.choices)
    return {}


def _workflow_ids() -> list[str]:
    # Import lazily to avoid a cli -> submission -> cli import cycle at startup.
    from geneset_extractors.cli import build_parser

    root = build_parser()
    workflows = _subparser_choices(root, "command").get("workflows")
    return sorted(_subparser_choices(workflows, "workflow_command")) if workflows else []


def _assay_type(identifier: str) -> str:
    for prefix, assay in (("rna", "rna_seq"), ("gtex", "rna_seq"), ("motrpac", "rna_seq"), ("atac", "atac_seq"), ("methylation", "methylation"), ("ptm", "proteomics"), ("proteomics", "proteomics"), ("splice", "splicing"), ("cnv", "copy_number"), ("drug", "drug_response"), ("morphology", "morphology"), ("hubmap", "anatomical_reference"), ("lincs", "perturbation"), ("calr", "calorimetry")):
        if identifier.startswith(prefix):
            return assay
    return "multi_assay"


def list_contracts() -> list[dict[str, Any]]:
    """Return a machine-readable view of current converter/workflow contracts."""
    contracts: list[dict[str, Any]] = []
    for identifier in sorted(CONVERTERS):
        try:
            spec = get_converter_spec(identifier)
            inputs = spec.get("inputs", [])
            outputs = spec.get("outputs", [])
            spec_error = None
        except Exception as exc:
            # Submission discovery must report every registered identifier even
            # when a pre-existing converter spec has an independent defect.
            inputs = []
            outputs = []
            spec_error = f"{type(exc).__name__}: {exc}"
        contracts.append(
            {
                "identifier": identifier,
                "kind": "converter",
                "source_module": get_converter(identifier).__name__,
                "cli_command": ["convert", identifier],
                "supported_assay_type": _assay_type(identifier),
                "expected_inputs": inputs,
                "expected_outputs": outputs,
                "test_fixture": "geneset_extractors/resources/submission_toy_deg.tsv" if identifier == "rna_deg" else None,
                "smoke_test_command": "geneset-extractors submission validate rna_deg" if identifier == "rna_deg" else None,
                "public_api_stability": "stable",
                "contract_metadata_error": spec_error,
            }
        )
    for identifier in _workflow_ids():
        contracts.append(
            {
                "identifier": identifier,
                "kind": "workflow",
                "source_module": f"geneset_extractors.workflows.{identifier}",
                "cli_command": ["workflows", identifier],
                "supported_assay_type": _assay_type(identifier),
                "expected_inputs": [{"name": "workflow-specific inputs", "type": "see CLI help"}],
                "expected_outputs": [{"name": "workflow artifacts", "type": "workflow-specific"}],
                "test_fixture": None,
                "smoke_test_command": None,
                "public_api_stability": "stable",
            }
        )
    return sorted(contracts, key=lambda item: (str(item["kind"]), str(item["identifier"])))


def describe_contract(identifier: str) -> dict[str, Any]:
    for contract in list_contracts():
        if contract["identifier"] == identifier:
            return contract
    raise KeyError(f"Unknown submission workflow/converter identifier: {identifier}")


def _cli_registration_exists(contract: dict[str, Any]) -> bool:
    if contract["kind"] == "converter":
        return str(contract["identifier"]) in CONVERTERS
    return str(contract["identifier"]) in set(_workflow_ids())


def _run_rna_deg_smoke() -> dict[str, Any]:
    from geneset_extractors.cli import main

    fixture = Path(__file__).parent / "resources" / "submission_toy_deg.tsv"
    if not fixture.is_file():
        return {"ok": False, "message": f"missing smoke fixture: {fixture}"}
    with tempfile.TemporaryDirectory(prefix="geneset_extractors_submission_") as temp_dir:
        out_dir = Path(temp_dir) / "rna_deg"
        code = main(
            [
                "convert", "rna_deg", "--deg_tsv", str(fixture), "--out_dir", str(out_dir),
                "--organism", "human", "--genome_build", "hg38", "--signature_name", "submission_smoke",
                "--top_k", "2", "--emit_gmt", "false", "--emit_full", "false",
            ]
        )
        required = ["geneset.tsv", "geneset.meta.json", "geneset.provenance.json"]
        missing = [name for name in required if not (out_dir / name).is_file()]
        if code != 0 or missing:
            return {"ok": False, "message": f"smoke converter failed code={code}; missing={missing}"}
        try:
            schema = Path(__file__).parent / "schemas" / "geneset_metadata.schema.json"
            validate_output_dir(out_dir, schema)
        except Exception as exc:
            return {"ok": False, "message": f"smoke output contract failed: {exc}"}
    return {"ok": True, "message": "rna_deg smoke fixture produced the standard final output contract"}


def validate_contract(identifier: str) -> dict[str, Any]:
    """Check registration/import and run an available, declared smoke test."""
    try:
        contract = describe_contract(identifier)
    except KeyError as exc:
        return {"identifier": identifier, "ok": False, "checks": {"registered": False}, "error": str(exc)}
    checks: dict[str, Any] = {"registered": _cli_registration_exists(contract)}
    try:
        importlib.import_module(str(contract["source_module"]))
    except Exception as exc:
        checks["module_importable"] = False
        checks["import_error"] = f"{type(exc).__name__}: {exc}"
    else:
        checks["module_importable"] = True
    smoke_command = contract.get("smoke_test_command")
    checks["declared_smoke_test"] = smoke_command is not None
    if identifier == "rna_deg" and checks["module_importable"]:
        smoke = _run_rna_deg_smoke()
        checks["smoke_test"] = smoke
    else:
        checks["smoke_test"] = {"ok": True, "message": "no declared low-cost smoke test"}
    ok = bool(checks["registered"] and checks["module_importable"] and checks["smoke_test"]["ok"])
    return {"identifier": identifier, "ok": ok, "contract": contract, "checks": checks}
