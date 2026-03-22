from __future__ import annotations

import json
from pathlib import Path


def load_provenance(out_dir: str | Path) -> dict[str, object]:
    return json.loads((Path(out_dir) / "geneset.provenance.json").read_text(encoding="utf-8"))


def file_node_for_role(provenance: dict[str, object], role: str) -> dict[str, object]:
    for node in provenance.get("nodes", []):
        if isinstance(node, dict) and node.get("kind") == "file" and node.get("role") == role:
            return node
    raise AssertionError(f"missing file provenance node for role={role}")


def assert_node_has_structured_resource_metadata(node: dict[str, object]) -> None:
    access = node.get("access", {})
    identifiers = node.get("identifiers", {})
    metadata = node.get("metadata", {})
    assert isinstance(access, dict)
    assert isinstance(identifiers, dict)
    assert isinstance(metadata, dict)
    assert access.get("local_path")
    assert (
        access.get("download_url")
        or access.get("landing_page_url")
        or access.get("canonical_uri")
        or identifiers.get("stable_id")
        or identifiers.get("persistent_id")
        or metadata.get("provider")
        or metadata.get("version")
        or metadata.get("license")
        or metadata.get("source")
    ), f"expected structured resource metadata beyond local_path for node={node}"


def assert_manifest_has_enriched_columns(rows: list[dict[str, str]]) -> None:
    assert rows
    required = {"path", "geneset_id", "label", "meta_path", "provenance_path", "focus_node_id"}
    assert required.issubset(rows[0].keys())
