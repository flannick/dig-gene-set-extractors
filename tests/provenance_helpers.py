from __future__ import annotations

import json
from pathlib import Path


def _graph_from_payload(payload: dict[str, object]) -> dict[str, object]:
    if "nodes" in payload and "edges" in payload:
        return payload
    if len(payload) != 1:
        raise AssertionError(f"expected exactly one provenance graph, found keys={list(payload.keys())}")
    graph = next(iter(payload.values()))
    if not isinstance(graph, dict):
        raise AssertionError("provenance graph payload is not an object")
    return graph


def load_provenance(out_dir: str | Path) -> dict[str, object]:
    payload = json.loads((Path(out_dir) / "geneset.provenance.json").read_text(encoding="utf-8"))
    return _graph_from_payload(payload)


def load_provenance_payload(out_dir: str | Path) -> dict[str, object]:
    return json.loads((Path(out_dir) / "geneset.provenance.json").read_text(encoding="utf-8"))


def file_node_for_role(provenance: dict[str, object], role: str) -> dict[str, object]:
    role_marker = f"(role: {role})"
    for node in provenance.get("nodes", []):
        if not isinstance(node, dict) or node.get("type") != "File":
            continue
        if role_marker in str(node.get("description", "")):
            return node
        if f":{role.lower()}:" in str(node.get("id", "")).lower():
            return node
    raise AssertionError(f"missing file provenance node for role={role}")


def assert_node_has_structured_resource_metadata(node: dict[str, object]) -> None:
    c2m2 = node.get("c2m2_properties", {})
    assert node.get("dcc_url")
    assert node.get("drc_url")
    assert isinstance(c2m2, dict)
    assert c2m2.get("local_id")
    assert (
        "resource" in str(node.get("description", "")).lower()
        or "provider=" in str(node.get("description", "")).lower()
        or "source=" in str(node.get("description", "")).lower()
        or str(node.get("dcc_url", "")).startswith(("http://", "https://", "gs://"))
        or str(node.get("drc_url", "")).startswith(("http://", "https://", "gs://"))
    ), f"expected structured resource metadata beyond local_path for node={node}"


def assert_manifest_has_enriched_columns(rows: list[dict[str, str]]) -> None:
    assert rows
    required = {"path", "geneset_id", "label", "meta_path", "provenance_path", "focus_node_id"}
    assert required.issubset(rows[0].keys())
