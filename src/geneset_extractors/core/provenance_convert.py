"""Convert persisted DIG legacy provenance graphs to DAPPER sidecars.

This module intentionally only loads and validates existing sidecars.  It does
not call metadata provenance builders, so conversion cannot rerun extraction
or replace the recorded legacy graph.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from geneset_extractors.core.dapper_provenance import build_dapper_provenance, write_dapper_provenance
from geneset_extractors.core.provenance import write_canonical_json


LEGACY_FILENAMES = ("geneset.provenance.legacy.json", "geneset.provenance.json")
DAPPER_FILENAME = "geneset.provenance.dapper.yaml"
_NODE_TYPES = {"File", "AnalysisType", "GeneSet", "GeneSetCollection"}
_EDGE_LABELS = {"data input", "metadata input", "data output"}


@dataclass(frozen=True)
class ConversionTarget:
    provenance: Path
    metadata: Path
    output: Path


def _load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"{label} is not valid JSON: {exc}") from exc
    except OSError as exc:
        raise ValueError(f"cannot read {label}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object")
    return value


def normalize_legacy_payload(payload: dict[str, Any], source: Path) -> dict[str, Any]:
    """Normalize the original single-graph DIG layout to the graph-map layout.

    Older DIG releases wrote one ``{nodes, edges}`` object directly; current
    releases wrap one or more graphs in a mapping keyed by gene-set id.
    Nothing else is reshaped, avoiding loss of provenance facts.
    """
    graphs = {"legacy": payload} if "nodes" in payload or "edges" in payload else payload
    if not graphs:
        raise ValueError(f"{source}: legacy provenance contains no graphs")
    for graph_key, graph in graphs.items():
        prefix = f"{source} graph {graph_key!r}"
        if not isinstance(graph, dict):
            raise ValueError(f"{prefix}: expected an object")
        nodes, edges = graph.get("nodes"), graph.get("edges")
        if not isinstance(nodes, list) or not nodes:
            raise ValueError(f"{prefix}: expected a non-empty nodes list")
        if not isinstance(edges, list):
            raise ValueError(f"{prefix}: expected an edges list")
        ids: set[str] = set()
        for index, node in enumerate(nodes):
            if not isinstance(node, dict):
                raise ValueError(f"{prefix}: node {index} is not an object")
            node_id, node_type = node.get("id"), node.get("type")
            if not isinstance(node_id, str) or not node_id:
                raise ValueError(f"{prefix}: node {index} has no non-empty id")
            if node_id in ids:
                raise ValueError(f"{prefix}: duplicate node id {node_id!r}")
            if node_type not in _NODE_TYPES:
                raise ValueError(f"{prefix}: unsupported node type {node_type!r} for {node_id!r}")
            ids.add(node_id)
        for index, edge in enumerate(edges):
            if not isinstance(edge, dict):
                raise ValueError(f"{prefix}: edge {index} is not an object")
            source_id, target_id, label = edge.get("source"), edge.get("target"), edge.get("label")
            if source_id not in ids or target_id not in ids:
                raise ValueError(f"{prefix}: edge {index} references a missing node")
            if label not in _EDGE_LABELS:
                raise ValueError(f"{prefix}: unsupported edge label {label!r}; cannot represent it in DAPPER")
    return graphs


def load_legacy_provenance(path: Path) -> dict[str, Any]:
    return normalize_legacy_payload(_load_json(path, f"legacy provenance {path}"), path)


def load_conversion_metadata(path: Path) -> dict[str, Any]:
    metadata = _load_json(path, f"metadata {path}")
    # These are the sections the existing mapper uses to classify the focus
    # node and preserve gene-set/GMT context. Older metadata may omit newer
    # optional fields, so do not require the current full metadata schema.
    for key in ("gene_set", "converter", "summary", "input", "output"):
        if not isinstance(metadata.get(key), dict):
            raise ValueError(f"metadata {path} is incompatible: expected object field {key!r}")
    return metadata


def discover_legacy_provenance(input_path: Path, recursive: bool) -> list[Path]:
    if input_path.is_file():
        return [input_path]
    if not input_path.is_dir():
        raise ValueError(f"input does not exist or is not a file/directory: {input_path}")
    candidates = (
        list(input_path.rglob("geneset.provenance*.json")) if recursive
        else [path for name in LEGACY_FILENAMES if (path := input_path / name).is_file()]
    )
    by_parent: dict[Path, Path] = {}
    for path in candidates:
        if path.name not in LEGACY_FILENAMES:
            continue
        prior = by_parent.get(path.parent)
        if prior is None or path.name == LEGACY_FILENAMES[0]:
            by_parent[path.parent] = path
    return sorted(by_parent.values())


def conversion_targets(
    input_path: Path, *, metadata: Path | None, output: Path | None, recursive: bool
) -> list[ConversionTarget]:
    provenance_paths = discover_legacy_provenance(input_path, recursive)
    if not provenance_paths:
        raise ValueError(f"no legacy provenance files found under {input_path}")
    if len(provenance_paths) > 1 and (metadata is not None or output is not None):
        raise ValueError("--metadata and --out require one input provenance file")
    return [
        ConversionTarget(
            provenance=provenance,
            metadata=metadata or provenance.parent / "geneset.meta.json",
            output=output or provenance.parent / DAPPER_FILENAME,
        )
        for provenance in provenance_paths
    ]


def convert_legacy_provenance(target: ConversionTarget, *, overwrite: bool) -> str:
    if target.output.exists() and not overwrite:
        return "skipped"
    if not target.metadata.is_file():
        raise ValueError(f"metadata file is missing: {target.metadata}")
    legacy = load_legacy_provenance(target.provenance)
    metadata = load_conversion_metadata(target.metadata)
    # The established writer supplies DAPPER-ID-1 minting, reference rewrites,
    # deterministic YAML, and optional GMT row export without editing the GMT.
    write_dapper_provenance(target.output, legacy, metadata)
    return "converted"


def duplicates_backup_path(path: Path) -> Path:
    return path.parent / "geneset.provenance.duplicates.json"


def deduplicate_legacy_provenance(path: Path, *, overwrite: bool) -> tuple[int, int]:
    """Back up and coalesce legacy File nodes identical under DIG's DAPPER map.

    The existing converter is deliberately the equivalence authority: nodes
    are coalesced only when it mints the same DAPPER class and ID. Every edge
    is redirected to the retained node, then exact source/target/label
    duplicates are collapsed. The original bytes remain in the companion
    ``geneset.provenance.duplicates.json`` backup.
    """
    payload = _load_json(path, f"legacy provenance {path}")
    graphs = normalize_legacy_payload(payload, path)
    metadata_path = path.parent / "geneset.meta.json"
    metadata = load_conversion_metadata(metadata_path) if metadata_path.is_file() else {}
    replacements_by_graph: dict[int, dict[str, str]] = {}
    removed_nodes = 0
    for graph in graphs.values():
        replacements: dict[str, str] = {}
        nodes = graph["nodes"]
        seen: dict[tuple[str, str], str] = {}
        retained: list[dict[str, Any]] = []
        for node in nodes:
            if node.get("type") != "File":
                retained.append(node)
                continue
            converted = build_dapper_provenance({"node": {"nodes": [node], "edges": []}}, metadata)
            bucket = "c2m2_files" if converted.get("c2m2_files") else "files"
            mapped_nodes = converted.get(bucket, [])
            if len(mapped_nodes) != 1 or not isinstance(mapped_nodes[0].get("id"), str):
                raise ValueError(f"{path}: could not map File node {node.get('id')!r} for deduplication")
            key = (bucket, mapped_nodes[0]["id"])
            node_id = str(node["id"])
            if key in seen:
                replacements[node_id] = seen[key]
                removed_nodes += 1
            else:
                seen[key] = node_id
                retained.append(node)
        graph["nodes"] = retained
        replacements_by_graph[id(graph)] = replacements
    if not any(replacements_by_graph.values()):
        return 0, 0
    removed_edges = 0
    for graph in graphs.values():
        replacements = replacements_by_graph[id(graph)]
        retained_edges: list[dict[str, Any]] = []
        seen_edges: set[tuple[object, object, object]] = set()
        for edge in graph["edges"]:
            rewritten = dict(edge)
            for field in ("source", "target"):
                if rewritten.get(field) in replacements:
                    rewritten[field] = replacements[rewritten[field]]
            key = (rewritten.get("source"), rewritten.get("target"), rewritten.get("label"))
            if key in seen_edges:
                removed_edges += 1
                continue
            seen_edges.add(key)
            retained_edges.append(rewritten)
        graph["edges"] = retained_edges
    backup = duplicates_backup_path(path)
    if backup.exists() and not overwrite:
        raise FileExistsError(f"duplicate provenance backup already exists: {backup}; pass --overwrite to replace it")
    backup.write_bytes(path.read_bytes())
    # Preserve the original layout (direct graph or graph map), only changing
    # the redundant nodes/edges required for DAPPER compatibility.
    write_canonical_json(path, graphs["legacy"] if "nodes" in payload or "edges" in payload else graphs)
    return removed_nodes, removed_edges
