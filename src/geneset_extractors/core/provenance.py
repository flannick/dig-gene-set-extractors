from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import sys
import threading
from typing import Any

from geneset_extractors.hashing import sha256_file, stable_hash_object


STANDARD_NAME = "dig.geneset"
STANDARD_VERSION = "0.1.0"
REPO_URL = "https://github.com/flannick/dig-gene-set-extractors"
_RUNTIME = threading.local()


@dataclass(frozen=True)
class RuntimeContext:
    converter_name: str
    entrypoint: str
    command: list[str]
    overlay_path: str | None


def activate_runtime_context(converter_name: str, overlay_path: str | None = None) -> None:
    command = list(sys.argv)
    entrypoint = f"geneset-extractors convert {converter_name}"
    if len(command) >= 3 and command[1] == "convert":
        entrypoint = " ".join(command[:3])
    _RUNTIME.context = RuntimeContext(
        converter_name=converter_name,
        entrypoint=entrypoint,
        command=command,
        overlay_path=(str(overlay_path).strip() if overlay_path else None),
    )


def get_runtime_context() -> RuntimeContext | None:
    ctx = getattr(_RUNTIME, "context", None)
    return ctx if isinstance(ctx, RuntimeContext) else None


def write_canonical_json(path: str | Path, payload: dict[str, Any]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def load_overlay(path: str | Path | None) -> dict[str, Any]:
    if not path:
        return {}
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("provenance overlay must be a JSON object")
    return payload


def stable_file_node_id(local_path: str, role: str, sha256: str | None) -> str:
    digest = stable_hash_object({"local_path": local_path, "role": role, "sha256": sha256 or ""})[:24]
    return f"file:{digest}"


def stable_operation_id(converter_name: str, geneset_id: str, input_ids: list[str]) -> str:
    digest = stable_hash_object(
        {
            "converter_name": converter_name,
            "geneset_id": geneset_id,
            "input_ids": sorted(input_ids),
        }
    )[:24]
    return f"operation:{digest}"


def geneset_node_id(geneset_id: str) -> str:
    return f"geneset:{geneset_id}"


def build_output_file_record(meta_dir: Path, file_record: dict[str, Any]) -> dict[str, Any]:
    raw_path = str(file_record.get("path", "")).strip()
    if not raw_path:
        raise ValueError("output file record missing path")
    path_obj = Path(raw_path)
    resolved = path_obj if path_obj.is_absolute() else (meta_dir / path_obj)
    sha256 = sha256_file(resolved) if resolved.exists() and resolved.is_file() else None
    size_bytes = resolved.stat().st_size if resolved.exists() and resolved.is_file() else None
    return {
        "path": raw_path,
        "role": str(file_record.get("role", "output_artifact")),
        "sha256": sha256,
        "size_bytes": size_bytes,
        "access_level": "local_only",
        "local_path": str(resolved),
    }


def _overlay_for_file(file_record: dict[str, Any], overlay: dict[str, Any]) -> dict[str, Any]:
    inputs = overlay.get("inputs", {})
    if not isinstance(inputs, dict):
        return {}
    role_key = f"role:{str(file_record.get('role', '')).strip()}"
    candidates = [role_key]
    raw_path = str(file_record.get("path", "")).strip()
    local_path = str(file_record.get("local_path", "")).strip()
    if raw_path:
        candidates.insert(0, raw_path)
    if local_path:
        candidates.insert(0, local_path)
        try:
            candidates.insert(0, str(Path(local_path).resolve()))
        except Exception:
            pass
    for key in candidates:
        value = inputs.get(key)
        if isinstance(value, dict):
            return value
    return {}


def merge_file_overlay(file_record: dict[str, Any], overlay: dict[str, Any]) -> dict[str, Any]:
    merged = dict(file_record)
    for key, value in _overlay_for_file(file_record, overlay).items():
        merged.setdefault(key, value)
    return merged


def build_file_node(file_record: dict[str, Any], overlay: dict[str, Any]) -> dict[str, Any]:
    merged = merge_file_overlay(file_record, overlay)
    local_path = str(merged.get("local_path") or merged.get("path") or "")
    role = str(merged.get("role", "input"))
    node_id = str(merged.get("node_id") or stable_file_node_id(local_path, role, merged.get("sha256")))
    hashes = {}
    if merged.get("sha256"):
        hashes["sha256"] = merged["sha256"]
    identifiers = {}
    if merged.get("persistent_id"):
        identifiers["persistent_id"] = merged["persistent_id"]
    if merged.get("stable_id"):
        identifiers["stable_id"] = merged["stable_id"]
    metadata = {}
    for key in ("provider", "version", "license", "source", "description"):
        value = merged.get(key)
        if value not in (None, ""):
            metadata[key] = value
    return {
        "id": node_id,
        "kind": "file",
        "label": str(merged.get("label") or Path(local_path or role).name or role),
        "role": role,
        "format": str(merged.get("format") or Path(local_path or merged.get("path", "")).suffix.lstrip(".") or "unknown"),
        "access": {
            "canonical_uri": merged.get("canonical_uri"),
            "download_url": merged.get("download_url"),
            "landing_page_url": merged.get("landing_page_url"),
            "local_path": local_path or None,
            "access_level": merged.get("access_level", "local_only"),
        },
        "identifiers": identifiers,
        "hashes": hashes,
        "size_bytes": merged.get("size_bytes"),
        "metadata": metadata,
        "extensions": merged.get("extensions", {}),
    }


def build_geneset_node(metadata_payload: dict[str, Any], overlay: dict[str, Any]) -> dict[str, Any]:
    gene_set = dict(metadata_payload.get("gene_set", {}))
    overlay_gene_set = overlay.get("gene_set", {})
    if isinstance(overlay_gene_set, dict):
        for key, value in overlay_gene_set.items():
            gene_set.setdefault(key, value)
    return {
        "id": str(gene_set["id"]),
        "kind": "geneset",
        "label": str(gene_set.get("name", metadata_payload.get("geneset_id", ""))),
        "assay": gene_set.get("assay"),
        "data_type": gene_set.get("data_type"),
        "organism": gene_set.get("organism"),
        "genome_build": gene_set.get("genome_build"),
        "n_genes": gene_set.get("n_genes"),
        "primary_artifact": gene_set.get("primary_artifact"),
        "metadata": {k: v for k, v in gene_set.items() if k not in {"id", "name", "assay", "data_type", "organism", "genome_build", "n_genes", "primary_artifact"}},
    }


def build_operation(
    metadata_payload: dict[str, Any],
    input_nodes: list[dict[str, Any]],
    output_nodes: list[dict[str, Any]],
    overlay: dict[str, Any],
) -> dict[str, Any]:
    converter = dict(metadata_payload.get("converter", {}))
    code = dict(converter.get("code", {}))
    execution = dict(converter.get("execution", {}))
    overlay_operation = overlay.get("operation", {})
    if isinstance(overlay_operation, dict):
        for key, value in overlay_operation.items():
            if key in {"script_url", "notebook_url"}:
                if code.get(key) in (None, ""):
                    code[key] = value
            elif key in {"container_image", "workspace_template_url"}:
                if execution.get(key) in (None, ""):
                    execution[key] = value
    focus_node_id = str(metadata_payload.get("provenance", {}).get("focus_node_id", ""))
    operation_id = stable_operation_id(
        str(converter.get("name", "unknown")),
        str(metadata_payload.get("geneset_id", "")),
        [str(node["id"]) for node in input_nodes],
    )
    command = execution.get("command")
    if not command:
        command = ["geneset-extractors", "convert", str(converter.get("name", "unknown")), "..."]
    return {
        "id": operation_id,
        "label": f"Extract {metadata_payload.get('gene_set', {}).get('name', metadata_payload.get('geneset_id', 'gene set'))}",
        "operation_type": "extract_gene_set",
        "method": str(converter.get("name", "unknown")),
        "inputs": [str(node["id"]) for node in input_nodes],
        "outputs": [focus_node_id] + [str(node["id"]) for node in output_nodes],
        "code": {
            "repo_url": code.get("repo_url", REPO_URL),
            "git_commit": code.get("git_commit"),
            "module": code.get("module"),
            "script_url": code.get("script_url"),
            "notebook_url": code.get("notebook_url"),
        },
        "parameters": converter.get("parameters", {}),
        "replay": {
            "mode": execution.get("mode", "cli"),
            "entrypoint": execution.get("entrypoint"),
            "command": command,
            "container_image": execution.get("container_image"),
            "workspace_template_url": execution.get("workspace_template_url"),
            "notes": execution.get("notes", "Template replay command; local inputs may require overlay URLs for portal replay."),
        },
    }


def build_edges(input_nodes: list[dict[str, Any]], operation_id: str, focus_node_id: str, output_nodes: list[dict[str, Any]]) -> list[dict[str, Any]]:
    edges: list[dict[str, Any]] = []
    for node in input_nodes:
        edges.append({"source": str(node["id"]), "target": operation_id, "kind": "used_by"})
    edges.append({"source": operation_id, "target": focus_node_id, "kind": "generated"})
    for node in output_nodes:
        edges.append({"source": focus_node_id, "target": str(node["id"]), "kind": "materialized_as"})
    return edges
