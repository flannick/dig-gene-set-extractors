from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
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
        {"path": "geneset.provenance.json", "role": "provenance_json"},
    ]
    for item in defaults:
        key = (item["path"], item["role"])
        if key not in seen:
            out.append(item)
    return out


def _build_provenance_payload(metadata_payload: dict[str, Any], meta_dir: Path, overlay_path: str | None) -> dict[str, Any]:
    runtime_ctx = get_runtime_context()
    overlay = load_overlay(overlay_path or (runtime_ctx.overlay_path if runtime_ctx else None))
    input_files = [dict(item) for item in metadata_payload.get("input", {}).get("files", [])]
    input_nodes = [build_file_node(item, overlay) for item in input_files]
    output_records = [build_output_file_record(meta_dir, item) for item in metadata_payload.get("output", {}).get("files", [])]
    output_nodes = [build_file_node(item, overlay) for item in output_records]
    geneset_node = build_geneset_node(metadata_payload, overlay)
    operation = build_operation(metadata_payload, input_nodes, output_nodes, overlay)
    return {
        "standard_name": STANDARD_NAME,
        "standard_version": STANDARD_VERSION,
        "file_type": "provenance",
        "focus_node_id": metadata_payload["provenance"]["focus_node_id"],
        "nodes": input_nodes + [geneset_node] + output_nodes,
        "operations": [operation],
        "edges": build_edges(input_nodes, str(operation["id"]), str(geneset_node["id"]), output_nodes),
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
