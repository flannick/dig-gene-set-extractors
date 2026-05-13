from __future__ import annotations

import base64
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import shlex
import sys
import threading
from typing import Any
from urllib.parse import unquote, urlparse
from uuid import NAMESPACE_URL, uuid5

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
    provenance_mirror_local_prefix: str | None
    provenance_mirror_remote_prefix: str | None


def _extract_cli_option(command: list[str], flag: str) -> str | None:
    for index, token in enumerate(command):
        if token == flag and index + 1 < len(command):
            value = str(command[index + 1]).strip()
            return value or None
        if token.startswith(flag + "="):
            value = token.split("=", 1)[1].strip()
            return value or None
    return None


def _normalize_cli_argv(argv: list[str]) -> list[str]:
    if not argv:
        return []
    tokens = [str(token) for token in argv]
    if len(tokens) >= 3 and tokens[1] == "-m" and tokens[2] == "geneset_extractors.cli":
        normalized = list(tokens)
        normalized[0] = "python"
        return normalized
    first = tokens[0]
    if first.endswith("/geneset_extractors/cli.py") or first.endswith("\\geneset_extractors\\cli.py") or first.endswith("/cli.py") or first.endswith("\\cli.py"):
        return ["python", "-m", "geneset_extractors.cli", *tokens[1:]]
    return tokens


def normalize_cli_command(command: list[str] | str | None) -> list[str] | str | None:
    if command is None:
        return None
    if isinstance(command, list):
        return _normalize_cli_argv(command)
    raw = str(command).strip()
    if not raw:
        return raw
    try:
        tokens = shlex.split(raw)
    except ValueError:
        return raw
    normalized = _normalize_cli_argv(tokens)
    return shlex.join(normalized)


def normalize_cli_entrypoint(entrypoint: str | None) -> str | None:
    if entrypoint in (None, ""):
        return entrypoint
    raw = str(entrypoint).strip()
    try:
        tokens = shlex.split(raw)
    except ValueError:
        tokens = raw.split()
    if not tokens:
        return raw
    if len(tokens) >= 3 and tokens[0] == "python" and tokens[1] == "-m" and tokens[2] == "geneset_extractors.cli":
        return " ".join(["geneset-extractors", *tokens[3:]]) or "geneset-extractors"
    if len(tokens) >= 3 and tokens[1] == "-m" and tokens[2] == "geneset_extractors.cli":
        return " ".join(["geneset-extractors", *tokens[3:]]) or "geneset-extractors"
    first = tokens[0]
    if first.endswith("/geneset_extractors/cli.py") or first.endswith("\\geneset_extractors\\cli.py") or first.endswith("/cli.py") or first.endswith("\\cli.py"):
        return " ".join(["geneset-extractors", *tokens[1:]]) or "geneset-extractors"
    return raw


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
        provenance_mirror_local_prefix=_extract_cli_option(command, "--provenance_mirror_local_prefix"),
        provenance_mirror_remote_prefix=_extract_cli_option(command, "--provenance_mirror_remote_prefix"),
    )


def get_runtime_context() -> RuntimeContext | None:
    ctx = getattr(_RUNTIME, "context", None)
    return ctx if isinstance(ctx, RuntimeContext) else None


def write_canonical_json(path: str | Path, payload: dict[str, Any]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


def load_overlay(path: str | Path | None) -> dict[str, Any]:
    if not path:
        return {}
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("provenance overlay must be a JSON object")
    return payload


def load_graph_payload(path: str | Path | None) -> dict[str, Any]:
    if not path:
        return {}
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("provenance graph payload must be a JSON object")
    return payload


def _slug(value: str) -> str:
    out = []
    for char in str(value):
        if char.isalnum():
            out.append(char.lower())
        else:
            out.append("_")
    slug = "".join(out).strip("_")
    while "__" in slug:
        slug = slug.replace("__", "_")
    return slug or "item"


def _stable_uuid(parts: dict[str, object]) -> str:
    return str(uuid5(NAMESPACE_URL, stable_hash_object(parts)))


def stable_file_node_id(local_path: str, role: str, sha256: str | None) -> str:
    file_uuid = _stable_uuid({"local_path": local_path, "role": role, "sha256": sha256 or ""})
    return f"file:{_slug(role)}:{file_uuid}"


def mirror_provenance_path(
    path: str | None,
    mirror_local_prefix: str | None,
    mirror_remote_prefix: str | None,
) -> str | None:
    raw_path = str(path or "").strip()
    if not raw_path or not mirror_local_prefix or not mirror_remote_prefix:
        return raw_path or None
    if "://" in raw_path and not raw_path.startswith("file://"):
        return raw_path
    try:
        local_root = Path(mirror_local_prefix).resolve()
        if raw_path.startswith("file://"):
            parsed = urlparse(raw_path)
            candidate = Path(unquote(parsed.path))
        else:
            candidate = Path(raw_path)
        resolved_candidate = candidate.resolve()
        relative = resolved_candidate.relative_to(local_root)
    except Exception:
        return raw_path
    remote_root = str(mirror_remote_prefix).rstrip("/")
    suffix = relative.as_posix()
    return f"{remote_root}/{suffix}" if suffix else remote_root


def _mirror_string_content(
    value: str,
    mirror_local_prefix: str | None,
    mirror_remote_prefix: str | None,
) -> str:
    if not value or not mirror_local_prefix or not mirror_remote_prefix:
        return value
    remote_root = str(mirror_remote_prefix).rstrip("/")
    local_root = str(Path(mirror_local_prefix).resolve()).rstrip("/")
    file_root = Path(local_root).as_uri().rstrip("/")
    for source_root in sorted({local_root, file_root}, key=len, reverse=True):
        value = value.replace(source_root, remote_root)
    return value


def _mirror_json_like(
    value: Any,
    mirror_local_prefix: str | None,
    mirror_remote_prefix: str | None,
) -> Any:
    if isinstance(value, str):
        return _mirror_string_content(value, mirror_local_prefix, mirror_remote_prefix)
    if isinstance(value, list):
        return [_mirror_json_like(item, mirror_local_prefix, mirror_remote_prefix) for item in value]
    if isinstance(value, dict):
        return {
            str(key): _mirror_json_like(item, mirror_local_prefix, mirror_remote_prefix)
            for key, item in value.items()
        }
    return value


def stable_operation_id(converter_name: str, geneset_id: str, input_ids: list[str]) -> str:
    op_uuid = _stable_uuid(
        {
            "converter_name": converter_name,
            "geneset_id": geneset_id,
            "input_ids": sorted(input_ids),
        }
    )
    return f"analysis:{_slug(converter_name)}:{op_uuid}"


def geneset_node_id(geneset_id: str) -> str:
    return f"geneset:{geneset_id}"


def _file_md5_base64(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return base64.b64encode(digest.digest()).decode("ascii")


def _normalize_md5(value: object, resolved: Path | None, *, local_path: str, role: str) -> str:
    if isinstance(value, str) and value.strip():
        raw = value.strip()
        try:
            digest = bytes.fromhex(raw)
        except ValueError:
            return raw
        return base64.b64encode(digest).decode("ascii")
    if resolved is not None and resolved.exists() and resolved.is_file():
        return _file_md5_base64(resolved)
    raise ValueError(f"Unable to compute MD5 for provenance file node role={role} path={local_path}")


def _path_to_uri(local_path: str) -> str | None:
    if not str(local_path).strip():
        return None
    try:
        return Path(local_path).resolve().as_uri()
    except Exception:
        return None


def _resolved_local_path_from_node(node: dict[str, Any]) -> Path | None:
    c2m2 = node.get("c2m2_properties", {})
    candidates: list[str] = []
    if isinstance(c2m2, dict):
        local_id = str(c2m2.get("local_id", "")).strip()
        if local_id:
            candidates.append(local_id)
    for key in ("dcc_url", "drc_url"):
        value = str(node.get(key, "")).strip()
        if value:
            candidates.append(value)
    for raw in candidates:
        try:
            if raw.startswith("file://"):
                parsed = urlparse(raw)
                candidate = Path(unquote(parsed.path))
            elif "://" in raw:
                continue
            else:
                candidate = Path(raw)
            resolved = candidate.resolve()
            if resolved.exists() and resolved.is_file():
                return resolved
        except Exception:
            continue
    return None


def _infer_node_description(label: str, role: str, kind: str, extra: str | None = None) -> str:
    base = f"{kind} node for {label} (role: {role})"
    if extra:
        return f"{base}. {extra}"
    return base


def _classify_input_edge(role: str) -> str:
    lowered = role.lower()
    metadata_tokens = (
        "gtf",
        "alias",
        "annotation",
        "metadata",
        "reference",
        "resource",
        "template",
        "manifest",
        "ontology",
    )
    if any(token in lowered for token in metadata_tokens):
        return "metadata input"
    return "data input"


def flatten_graph_payload(payload: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    nodes: list[dict[str, Any]] = []
    edges: list[dict[str, Any]] = []
    for graph in payload.values():
        if not isinstance(graph, dict):
            continue
        raw_nodes = graph.get("nodes", [])
        raw_edges = graph.get("edges", [])
        if isinstance(raw_nodes, list):
            nodes.extend(node for node in raw_nodes if isinstance(node, dict))
        if isinstance(raw_edges, list):
            edges.extend(edge for edge in raw_edges if isinstance(edge, dict))
    return nodes, edges


def merge_graph_components(
    node_groups: list[list[dict[str, Any]]],
    edge_groups: list[list[dict[str, Any]]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    merged_nodes: list[dict[str, Any]] = []
    seen_node_ids: set[str] = set()
    for group in node_groups:
        for node in group:
            node_id = str(node.get("id", "")).strip()
            if not node_id or node_id in seen_node_ids:
                continue
            seen_node_ids.add(node_id)
            merged_nodes.append(node)

    merged_edges: list[dict[str, Any]] = []
    seen_edge_ids: set[str] = set()
    for group in edge_groups:
        for edge in group:
            edge_id = str(edge.get("id", "")).strip()
            if not edge_id or edge_id in seen_edge_ids:
                continue
            seen_edge_ids.add(edge_id)
            merged_edges.append(edge)
    return merged_nodes, merged_edges


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
        "md5": file_record.get("md5"),
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
def build_file_node(
    file_record: dict[str, Any],
    overlay: dict[str, Any],
    *,
    mirror_local_prefix: str | None = None,
    mirror_remote_prefix: str | None = None,
) -> dict[str, Any]:
    merged = merge_file_overlay(file_record, overlay)
    local_path = str(merged.get("local_path") or merged.get("path") or "")
    resolved = Path(local_path).resolve() if local_path else None
    role = str(merged.get("role", "input"))
    sha256 = merged.get("sha256")
    mirrored_local_id = mirror_provenance_path(local_path, mirror_local_prefix, mirror_remote_prefix)
    mirrored_explicit_local_id = mirror_provenance_path(
        str(merged.get("local_id", "")).strip() or None,
        mirror_local_prefix,
        mirror_remote_prefix,
    )
    identity_path = (
        mirrored_explicit_local_id
        or str(merged.get("canonical_uri") or "").strip()
        or mirrored_local_id
        or local_path
    )
    node_id = str(merged.get("node_id") or stable_file_node_id(identity_path, role, sha256 if isinstance(sha256, str) else None))
    file_uuid = str(merged.get("_uuid") or _stable_uuid({"node_id": node_id, "role": role}))
    persistent_id = str(merged.get("persistent_id") or _stable_uuid({"node_id": node_id, "kind": "persistent"}))
    filename = str(
        merged.get("filename")
        or (Path(local_path).name if local_path else Path(str(merged.get("path", role))).name)
        or role
    )
    local_id = str(
        mirrored_explicit_local_id
        or merged.get("canonical_uri")
        or mirrored_local_id
        or merged.get("path")
        or filename
    )
    size_bytes = merged.get("size_bytes")
    if size_bytes is None and resolved is not None and resolved.exists() and resolved.is_file():
        size_bytes = resolved.stat().st_size
    md5 = _normalize_md5(
        merged.get("md5"),
        resolved if resolved is not None else None,
        local_path=local_path,
        role=role,
    )
    dcc_url = (
        merged.get("dcc_url")
        or merged.get("landing_page_url")
        or merged.get("download_url")
        or merged.get("canonical_uri")
        or mirrored_local_id
        or _path_to_uri(local_path)
        or f"urn:file:{file_uuid}"
    )
    drc_url = (
        merged.get("drc_url")
        or merged.get("canonical_uri")
        or merged.get("download_url")
        or merged.get("landing_page_url")
        or mirrored_local_id
        or _path_to_uri(local_path)
        or f"urn:file:{persistent_id}"
    )
    description_bits = []
    for key in ("provider", "version", "license", "source"):
        value = merged.get(key)
        if value not in (None, ""):
            description_bits.append(f"{key}={value}")
    return {
        "id": node_id,
        "type": "File",
        "name": str(merged.get("label") or filename),
        "description": _infer_node_description(
            str(merged.get("label") or filename),
            role,
            "File",
            "; ".join(description_bits) if description_bits else None,
        ),
        "dcc_url": str(dcc_url),
        "drc_url": str(drc_url),
        "c2m2_properties": {
            "filename": filename,
            "persistent_id": persistent_id,
            "local_id": local_id,
            "size_in_bytes": int(size_bytes or 0),
            "_uuid": file_uuid,
            "md5": md5,
        },
    }


def build_analysis_node(
    *,
    analysis_id: str,
    method: str,
    name: str,
    description: str,
    parameters: dict[str, Any],
    command: list[str] | str | None,
    entrypoint: str | None,
    observed_command: list[str] | str | None = None,
    repo_url: str | None = None,
    module: str | None = None,
    script_url: str | None = None,
    version: str | None = None,
    container_image: str | None = None,
    workspace_template_url: str | None = None,
    dcc_url: str | None = None,
    drc_url: str | None = None,
) -> dict[str, Any]:
    normalized_command = normalize_cli_command(command)
    normalized_observed_command = normalize_cli_command(observed_command)
    normalized_entrypoint = normalize_cli_entrypoint(entrypoint)
    if isinstance(normalized_command, list):
        command_text = shlex.join([str(token) for token in normalized_command])
    elif normalized_command in (None, ""):
        command_text = None
    else:
        command_text = str(normalized_command)
    resolved_repo_url = repo_url or REPO_URL
    resolved_script_url = script_url or resolved_repo_url
    return {
        "id": analysis_id,
        "type": "AnalysisType",
        "name": name,
        "description": description,
        "dcc_url": str(dcc_url or resolved_repo_url),
        "drc_url": str(drc_url or resolved_repo_url),
        "c2m2_properties": {
            "id": analysis_id,
            "name": name,
            "description": description,
            "synonyms": [],
        },
        "analysis": {
            "script_url": str(resolved_script_url),
            "version": version,
            "command": command_text,
            "observed_command": (
                shlex.join([str(token) for token in normalized_observed_command])
                if isinstance(normalized_observed_command, list)
                else (str(normalized_observed_command) if normalized_observed_command not in (None, "") else None)
            ),
            "parameters": parameters,
            "environment": {
                "mode": "cli",
                "entrypoint": normalized_entrypoint,
                "container_image": container_image,
                "workspace_template_url": workspace_template_url,
                "repo_url": resolved_repo_url,
                "module": module,
            },
        },
    }


def mirror_graph_payload(
    payload: dict[str, Any],
    mirror_local_prefix: str | None,
    mirror_remote_prefix: str | None,
) -> dict[str, Any]:
    if not mirror_local_prefix or not mirror_remote_prefix:
        return payload
    rewritten: dict[str, Any] = {}
    for graph_key, graph in payload.items():
        if not isinstance(graph, dict):
            rewritten[str(graph_key)] = graph
            continue
        id_map: dict[str, str] = {}
        new_nodes: list[dict[str, Any]] = []
        for node in graph.get("nodes", []):
            if not isinstance(node, dict):
                continue
            node_copy = json.loads(json.dumps(node))
            if node_copy.get("type") == "File":
                c2m2 = node_copy.get("c2m2_properties", {})
                if not isinstance(c2m2, dict):
                    c2m2 = {}
                    node_copy["c2m2_properties"] = c2m2
                role = "file"
                description = str(node_copy.get("description", ""))
                marker = "(role:"
                if marker in description:
                    role = description.split(marker, 1)[1].split(")", 1)[0].strip() or role
                original_identity = (
                    str(c2m2.get("local_id") or "").strip()
                    or str(node_copy.get("dcc_url") or "").strip()
                    or str(node_copy.get("drc_url") or "").strip()
                )
                mirrored_identity = mirror_provenance_path(
                    original_identity,
                    mirror_local_prefix,
                    mirror_remote_prefix,
                )
                if mirrored_identity and mirrored_identity != original_identity:
                    original_id = str(node_copy.get("id", ""))
                    sha256 = None
                    resolved_local = _resolved_local_path_from_node(node)
                    if resolved_local is not None:
                        try:
                            sha256 = sha256_file(resolved_local)
                        except Exception:
                            sha256 = None
                    new_id = stable_file_node_id(mirrored_identity, role, sha256)
                    id_map[original_id] = new_id
                    node_copy["id"] = new_id
                    file_uuid = _stable_uuid({"node_id": new_id, "role": role})
                    c2m2["_uuid"] = file_uuid
                    c2m2["persistent_id"] = _stable_uuid({"node_id": new_id, "kind": "persistent"})
                    c2m2["local_id"] = mirrored_identity
                    node_copy["dcc_url"] = _mirror_string_content(str(node_copy.get("dcc_url", "")), mirror_local_prefix, mirror_remote_prefix)
                    node_copy["drc_url"] = _mirror_string_content(str(node_copy.get("drc_url", "")), mirror_local_prefix, mirror_remote_prefix)
            elif node_copy.get("type") == "AnalysisType":
                analysis = node_copy.get("analysis", {})
                if isinstance(analysis, dict):
                    if isinstance(analysis.get("command"), str):
                        mirrored_command = _mirror_string_content(
                            str(analysis["command"]),
                            mirror_local_prefix,
                            mirror_remote_prefix,
                        )
                        analysis["command"] = normalize_cli_command(mirrored_command)
                    if isinstance(analysis.get("observed_command"), str):
                        mirrored_observed = _mirror_string_content(
                            str(analysis["observed_command"]),
                            mirror_local_prefix,
                            mirror_remote_prefix,
                        )
                        analysis["observed_command"] = normalize_cli_command(mirrored_observed)
                    if "parameters" in analysis:
                        analysis["parameters"] = _mirror_json_like(
                            analysis.get("parameters"),
                            mirror_local_prefix,
                            mirror_remote_prefix,
                        )
                    environment = analysis.get("environment", {})
                    if isinstance(environment, dict):
                        for key in ("entrypoint", "workspace_template_url"):
                            if isinstance(environment.get(key), str):
                                environment[key] = _mirror_string_content(
                                    str(environment[key]),
                                    mirror_local_prefix,
                                    mirror_remote_prefix,
                                )
                        if isinstance(environment.get("entrypoint"), str):
                            environment["entrypoint"] = normalize_cli_entrypoint(environment.get("entrypoint"))
            new_nodes.append(node_copy)
        new_edges: list[dict[str, Any]] = []
        for edge in graph.get("edges", []):
            if not isinstance(edge, dict):
                continue
            edge_copy = json.loads(json.dumps(edge))
            source = str(edge_copy.get("source", ""))
            target = str(edge_copy.get("target", ""))
            if source in id_map:
                edge_copy["source"] = id_map[source]
            if target in id_map:
                edge_copy["target"] = id_map[target]
            edge_copy["id"] = _mirror_string_content(str(edge_copy.get("id", "")), mirror_local_prefix, mirror_remote_prefix)
            new_edges.append(edge_copy)
        rewritten[str(graph_key)] = {
            "nodes": new_nodes,
            "edges": new_edges,
        }
    return rewritten


def build_geneset_node(metadata_payload: dict[str, Any], overlay: dict[str, Any]) -> dict[str, Any]:
    gene_set = dict(metadata_payload.get("gene_set", {}))
    overlay_gene_set = overlay.get("gene_set", {})
    if isinstance(overlay_gene_set, dict):
        for key, value in overlay_gene_set.items():
            gene_set.setdefault(key, value)
    geneset_id = str(gene_set["id"])
    name = str(gene_set.get("name", metadata_payload.get("geneset_id", geneset_id)))
    description = str(
        gene_set.get("description")
        or f"Derived gene set for converter={metadata_payload.get('converter', {}).get('name', 'unknown')}"
    )
    dcc_url = str(gene_set.get("dcc_url") or f"urn:geneset:{metadata_payload.get('geneset_id', geneset_id)}")
    drc_url = str(gene_set.get("drc_url") or dcc_url)
    return {
        "id": geneset_id,
        "type": "GeneSet",
        "name": name,
        "description": description,
        "dcc_url": dcc_url,
        "drc_url": drc_url,
        "c2m2_properties": {
            "id": geneset_id,
            "name": name,
            "description": description,
        },
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
            elif key in {"dcc_url", "drc_url", "description", "synonyms"}:
                code.setdefault(key, value)
    operation_id = stable_operation_id(
        str(converter.get("name", "unknown")),
        str(metadata_payload.get("geneset_id", "")),
        [str(node["id"]) for node in input_nodes],
    )
    method = str(converter.get("name", "unknown"))
    label = f"generate_{_slug(metadata_payload.get('geneset_id', method))}"
    analysis_description = str(
        code.get("description")
        or f"Analysis step that derives the gene set using converter {method}."
    )
    operation = build_analysis_node(
        analysis_id=operation_id,
        method=method,
        name=label,
        description=analysis_description,
        parameters=converter.get("parameters", {}),
        command=execution.get("command"),
        entrypoint=execution.get("entrypoint"),
        observed_command=execution.get("observed_command"),
        repo_url=code.get("repo_url", REPO_URL),
        module=code.get("module"),
        script_url=code.get("script_url") or code.get("notebook_url") or REPO_URL,
        version=code.get("git_commit"),
        container_image=execution.get("container_image"),
        workspace_template_url=execution.get("workspace_template_url"),
        dcc_url=code.get("dcc_url") or REPO_URL,
        drc_url=code.get("drc_url") or REPO_URL,
    )
    operation["c2m2_properties"]["synonyms"] = list(code.get("synonyms", [])) if isinstance(code.get("synonyms"), list) else []
    operation["analysis"]["environment"]["n_input_nodes"] = len(input_nodes)
    operation["analysis"]["environment"]["n_output_nodes"] = len(output_nodes)
    return operation


def build_edges(
    input_nodes: list[dict[str, Any]],
    operation_id: str,
    focus_node_id: str,
    output_nodes: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    edges: list[dict[str, Any]] = []
    for node in input_nodes:
        role = str(node.get("description", ""))
        label = _classify_input_edge(role)
        edge_id = f"{node['id']}_to_{operation_id}"
        edges.append(
            {
                "id": edge_id,
                "source": str(node["id"]),
                "target": operation_id,
                "label": label,
                "description": f"{node['name']} is a {label} to {operation_id}.",
            }
        )
    edges.append(
        {
            "id": f"{operation_id}_to_{focus_node_id}",
            "source": operation_id,
            "target": focus_node_id,
            "label": "data output",
            "description": f"{operation_id} produces {focus_node_id}.",
        }
    )
    for node in output_nodes:
        edges.append(
            {
                "id": f"{operation_id}_to_{node['id']}",
                "source": operation_id,
                "target": str(node["id"]),
                "label": "data output",
                "description": f"{operation_id} materializes output artifact {node['name']}.",
            }
        )
    return edges
