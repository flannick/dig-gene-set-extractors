from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import json
import mimetypes
from pathlib import Path
import shlex
import subprocess
import sys
import threading
from typing import Any

from geneset_extractors import __version__
from geneset_extractors.hashing import sha256_file, stable_hash_object
from geneset_extractors.core.provenance import (
    REPO_URL,
    STANDARD_NAME,
    STANDARD_VERSION,
    build_edges,
    build_file_node,
    build_geneset_node,
    build_operation,
    build_output_file_record,
    geneset_node_id,
    get_runtime_context,
    load_overlay,
    write_canonical_json,
)


def build_geneset_id(converter_name: str, file_hashes: list[str], params: dict[str, object]) -> str:
    return stable_hash_object(
        {
            "converter": converter_name,
            "file_hashes": sorted(file_hashes),
            "params_hash": stable_hash_object(params),
        }
    )[:24]


def write_metadata(path: str | Path, payload: dict[str, object]) -> None:
    p = Path(path)
    clean_payload = {k: v for k, v in payload.items() if not str(k).startswith("_")}
    write_canonical_json(p, clean_payload)
    if p.name == "geneset.meta.json":
        overlay_path = payload.get("_provenance_overlay_json")
        provenance_payload = _build_provenance_payload(clean_payload, p.parent, overlay_path if isinstance(overlay_path, str) else None)
        write_canonical_json(p.parent / "geneset.provenance.json", provenance_payload)


_GIT_COMMIT_CACHE: str | None = None
_GIT_COMMIT_LOCK = threading.Lock()
_INVOCATION_CONTEXT = threading.local()


def _resolve_git_commit() -> str:
    global _GIT_COMMIT_CACHE
    with _GIT_COMMIT_LOCK:
        if _GIT_COMMIT_CACHE is not None:
            return _GIT_COMMIT_CACHE
        repo_root = Path(__file__).resolve().parents[3]
        commit = "unknown"
        try:
            proc = subprocess.run(
                ["git", "rev-parse", "HEAD"],
                cwd=str(repo_root),
                check=False,
                capture_output=True,
                text=True,
            )
            if proc.returncode == 0:
                candidate = proc.stdout.strip()
                if candidate:
                    commit = candidate
        except Exception:
            commit = "unknown"
        _GIT_COMMIT_CACHE = commit
        return commit


def _default_gene_set_name(converter_name: str, parameters: dict[str, object], geneset_id: str) -> str:
    signature = str(parameters.get("signature_name", "")).strip()
    comparison = str(parameters.get("comparison_label", "")).strip()
    dataset = str(parameters.get("dataset_label", "")).strip()
    program = str(parameters.get("program_id") or parameters.get("program") or "").strip()
    sample_id = str(parameters.get("sample_id", "")).strip()
    pieces = [part for part in [dataset, signature or program, comparison or sample_id] if part]
    if pieces:
        return " | ".join(pieces)
    if signature:
        return signature
    return f"{converter_name}:{geneset_id}"


def _default_focus_node_id(geneset_id: str) -> str:
    return geneset_node_id(geneset_id)


def _ensure_output_files(output_files: list[dict[str, object]] | None) -> list[dict[str, object]]:
    out = [dict(item) for item in (output_files or [])]
    seen = {(str(item.get("path", "")), str(item.get("role", ""))) for item in out}
    defaults = [
        {"path": "geneset.tsv", "role": "selected_program"},
        {"path": "geneset.meta.json", "role": "metadata_json"},
    ]
    for item in defaults:
        key = (item["path"], item["role"])
        if key not in seen:
            out.append(item)
    return out


def _dedupe_output_files(meta_dir: Path, output_files: list[dict[str, object]]) -> list[dict[str, object]]:
    deduped: list[dict[str, object]] = []
    seen: set[tuple[str, str]] = set()
    for item in output_files:
        raw_path = str(item.get("path", "")).strip()
        role = str(item.get("role", "")).strip()
        if not raw_path:
            continue
        path_obj = Path(raw_path)
        resolved = path_obj if path_obj.is_absolute() else (meta_dir / path_obj)
        key = (str(resolved.resolve()), role)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(dict(item))
    return deduped


def _build_provenance_payload(metadata_payload: dict[str, Any], meta_dir: Path, overlay_path: str | None) -> dict[str, Any]:
    runtime_ctx = get_runtime_context()
    overlay = load_overlay(overlay_path or (runtime_ctx.overlay_path if runtime_ctx else None))
    input_files = [dict(item) for item in metadata_payload.get("input", {}).get("files", [])]
    input_nodes = [build_file_node(item, overlay) for item in input_files]
    output_files = _dedupe_output_files(meta_dir, [dict(item) for item in metadata_payload.get("output", {}).get("files", [])])
    output_records = [build_output_file_record(meta_dir, item) for item in output_files]
    output_nodes = [build_file_node(item, overlay) for item in output_records]
    geneset_node = build_geneset_node(metadata_payload, overlay)
    operation = build_operation(metadata_payload, input_nodes, output_nodes, overlay)
    graph_key = str(metadata_payload.get("geneset_id", metadata_payload["provenance"]["focus_node_id"]))
    return {
        graph_key: {
            "nodes": input_nodes + [geneset_node, operation] + output_nodes,
            "edges": build_edges(input_nodes, str(operation["id"]), str(geneset_node["id"]), output_nodes),
        }
    }


@contextmanager
def invocation_context(command_argv: list[str], cwd: str | Path | None = None):
    previous = getattr(_INVOCATION_CONTEXT, "value", None)
    _INVOCATION_CONTEXT.value = {
        "argv": [str(token) for token in command_argv],
        "cwd": str(cwd) if cwd is not None else str(Path.cwd()),
        "kind": "captured_cli_argv",
    }
    try:
        yield
    finally:
        _INVOCATION_CONTEXT.value = previous


def _current_invocation_context() -> dict[str, object] | None:
    value = getattr(_INVOCATION_CONTEXT, "value", None)
    return dict(value) if isinstance(value, dict) else None


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


def _command_flag(name: str) -> str:
    return f"--{name}"


def _append_parameter_tokens(tokens: list[str], name: str, value: object) -> None:
    flag = _command_flag(name)
    if value is None:
        return
    if isinstance(value, bool):
        tokens.extend([flag, "true" if value else "false"])
        return
    if isinstance(value, dict):
        tokens.extend([flag, json.dumps(value, sort_keys=True)])
        return
    if isinstance(value, (list, tuple)):
        if not value:
            return
        for item in value:
            _append_parameter_tokens(tokens, name, item)
        return
    tokens.extend([flag, str(value)])


def _reconstruct_command_argv(converter_name: str, parameters: dict[str, object]) -> list[str]:
    argv = [sys.executable, "-m", "geneset_extractors.cli", "convert", converter_name]
    for key in sorted(parameters):
        _append_parameter_tokens(argv, key, parameters[key])
    return argv


def _guess_media_type(path: str) -> str | None:
    media_type, _ = mimetypes.guess_type(path)
    return media_type


def _build_lineage_file_node(
    record: dict[str, object],
    direction: str,
    index: int,
) -> dict[str, object]:
    path = str(record.get("path", ""))
    role = str(record.get("role", "")) or f"{direction}_{index}"
    node: dict[str, object] = {
        "id": f"file:{direction}:{_slug(role)}:{index}",
        "node_type": "file",
        "direction": direction,
        "role": role,
        "path": path,
    }
    sha256 = record.get("sha256")
    if isinstance(sha256, str) and sha256:
        node["sha256"] = sha256
    else:
        candidate = Path(path)
        if path and candidate.exists() and candidate.is_file():
            node["sha256"] = sha256_file(candidate)
    media_type = _guess_media_type(path)
    if media_type:
        node["media_type"] = media_type
    return node


def _build_lineage(
    converter_name: str,
    parameters: dict[str, object],
    files: list[dict[str, str]],
    output_files: list[dict[str, str]] | None,
) -> dict[str, object]:
    input_nodes = [_build_lineage_file_node(record, "input", idx) for idx, record in enumerate(files, start=1)]
    output_nodes = [
        _build_lineage_file_node(record, "output", idx)
        for idx, record in enumerate(output_files or [], start=1)
    ]

    invocation = _current_invocation_context()
    if invocation is None:
        command_argv = _reconstruct_command_argv(converter_name, parameters)
        command_kind = "reconstructed_from_parameters"
        working_directory = str(Path.cwd())
    else:
        command_argv = [str(token) for token in invocation.get("argv", [])]
        command_kind = str(invocation.get("kind", "captured_cli_argv"))
        working_directory = str(invocation.get("cwd", Path.cwd()))

    process_id = "process:converter_invocation"
    processes = [
        {
            "id": process_id,
            "name": f"convert {converter_name}",
            "process_type": "converter_invocation",
            "command_kind": command_kind,
            "command_argv": command_argv,
            "command": shlex.join(command_argv),
            "working_directory": working_directory,
            "entrypoint": {"kind": "python_module", "module": "geneset_extractors.cli"},
            "parameters": parameters,
            "code": {"git_commit": _resolve_git_commit()},
        }
    ]

    edges = []
    if output_nodes:
        for input_node in input_nodes:
            for output_node in output_nodes:
                pair_hash = stable_hash_object(
                    {
                        "process_id": process_id,
                        "source": input_node["id"],
                        "target": output_node["id"],
                    }
                )[:12]
                edges.append(
                    {
                        "id": f"edge:{pair_hash}",
                        "source": input_node["id"],
                        "target": output_node["id"],
                        "process_id": process_id,
                        "edge_type": "process_step",
                    }
                )

    return {
        "graph_version": "1.0.0",
        "nodes": input_nodes + output_nodes,
        "edges": edges,
        "processes": processes,
    }


def make_metadata(
    converter_name: str,
    parameters: dict[str, object],
    data_type: str,
    assay: str,
    organism: str,
    genome_build: str,
    files: list[dict[str, str]],
    gene_annotation: dict[str, object],
    weights: dict[str, object],
    summary: dict[str, object],
    program_extraction: dict[str, object] | None = None,
    output_files: list[dict[str, str]] | None = None,
    gmt: dict[str, object] | None = None,
    provenance_overlay_json: str | None = None,
) -> dict[str, object]:
    file_hashes = [f["sha256"] for f in files]
    geneset_id = build_geneset_id(converter_name, file_hashes, parameters)
    focus_node_id = _default_focus_node_id(geneset_id)
    runtime_ctx = get_runtime_context()
    output_files_payload = _ensure_output_files(output_files)
    gene_set_name = _default_gene_set_name(converter_name, parameters, geneset_id)
    payload: dict[str, object] = {
        "standard_name": STANDARD_NAME,
        "standard_version": STANDARD_VERSION,
        "file_type": "metadata",
        "schema_version": "1.0.0",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "geneset_id": geneset_id,
        "gene_set": {
            "id": focus_node_id,
            "name": gene_set_name,
            "assay": assay,
            "data_type": data_type,
            "organism": organism,
            "genome_build": genome_build,
            "n_genes": int(summary.get("n_genes", 0)),
            "primary_artifact": {
                "path": "geneset.tsv",
                "role": "selected_program",
            },
        },
        "converter": {
            "name": converter_name,
            "version": __version__,
            "parameters": parameters,
            "code": {
                "repo_url": REPO_URL,
                "git_commit": _resolve_git_commit(),
                "module": f"geneset_extractors.extractors.converters.{converter_name}",
                "script_url": None,
                "notebook_url": None,
            },
            "execution": {
                "mode": "cli",
                "entrypoint": (
                    runtime_ctx.entrypoint if runtime_ctx is not None else f"geneset-extractors convert {converter_name}"
                ),
                "command": (
                    runtime_ctx.command
                    if runtime_ctx is not None and runtime_ctx.command
                    else ["geneset-extractors", "convert", converter_name, "..."]
                ),
                "container_image": None,
                "workspace_template_url": None,
                "notes": "Collapsed extraction operation for one emitted gene set.",
            },
        },
        "input": {
            "data_type": data_type,
            "assay": assay,
            "organism": organism,
            "genome_build": genome_build,
            "files": files,
        },
        "gene_annotation": gene_annotation,
        "weights": weights,
        "summary": summary,
        "provenance": {
            "path": "geneset.provenance.json",
            "focus_node_id": focus_node_id,
        },
    }
    if program_extraction is not None:
        payload["program_extraction"] = program_extraction
    payload["output"] = {"files": output_files_payload}
    if gmt is not None:
        payload["gmt"] = gmt
    payload["_provenance_overlay_json"] = provenance_overlay_json
    payload["lineage"] = _build_lineage(
        converter_name=converter_name,
        parameters=parameters,
        files=files,
        output_files=output_files_payload,
    )
    return payload


def input_file_record(
    path: str | Path,
    role: str,
    *,
    resource_record: dict[str, object] | None = None,
    **overrides: object,
) -> dict[str, object]:
    p = Path(path)
    payload: dict[str, object] = {
        "path": str(p),
        "local_path": str(p),
        "role": role,
        "sha256": sha256_file(p),
        "size_bytes": p.stat().st_size,
        "access_level": "local_only",
    }
    if resource_record is not None:
        payload.update(
            {
                "canonical_uri": resource_record.get("canonical_uri") or resource_record.get("url"),
                "download_url": resource_record.get("download_url") or resource_record.get("url"),
                "landing_page_url": resource_record.get("landing_page_url"),
                "persistent_id": resource_record.get("persistent_id"),
                "provider": resource_record.get("provider"),
                "version": resource_record.get("version"),
                "license": resource_record.get("license"),
                "stable_id": resource_record.get("stable_id") or resource_record.get("id"),
                "source": "resource_manager",
            }
        )
        if resource_record.get("url"):
            payload["access_level"] = "public"
        if resource_record.get("sha256"):
            payload["sha256"] = resource_record["sha256"]
    for key, value in overrides.items():
        if value is not None:
            payload[key] = value
    return payload


def enrich_manifest_row(out_dir: str | Path, group_dir: str | Path, row: dict[str, object]) -> dict[str, object]:
    out_root = Path(out_dir)
    group_path = Path(group_dir)
    meta_path = group_path / "geneset.meta.json"
    payload = json.loads(meta_path.read_text(encoding="utf-8"))
    provenance = payload.get("provenance", {})
    enriched = dict(row)
    enriched.setdefault("geneset_id", payload.get("geneset_id", ""))
    enriched.setdefault("label", payload.get("gene_set", {}).get("name", ""))
    enriched.setdefault("path", str(group_path.relative_to(out_root)))
    enriched.setdefault("meta_path", str(meta_path.relative_to(out_root)))
    enriched.setdefault(
        "provenance_path",
        str((group_path / str(provenance.get("path", "geneset.provenance.json"))).relative_to(out_root)),
    )
    enriched.setdefault("focus_node_id", provenance.get("focus_node_id", ""))
    return enriched
