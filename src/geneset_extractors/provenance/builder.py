from __future__ import annotations

import json
import re
import shlex
import sys
from pathlib import Path
from typing import Any

from geneset_extractors.hashing import stable_hash_object

from .checksums import checksum_list, file_size
from .environment import git_info, pip_freeze_text, python_environment, r_package_version, r_session_info_text
from .overlay import mirror_for_role, software_overlay, source_for_role


ROLE_RE = re.compile(r"\(role:\s*([^)]+)\)")

CONTEXT = {
    "schema": "https://schema.org/",
    "prov": "http://www.w3.org/ns/prov#",
    "dig": "https://dig.example.org/provenance#",
}


def _role_from_compact_node(node: dict[str, Any]) -> str:
    description = str(node.get("description", ""))
    match = ROLE_RE.search(description)
    if match:
        return match.group(1).strip()
    node_id = str(node.get("id", ""))
    parts = node_id.split(":")
    if len(parts) >= 2:
        return parts[1]
    return "file"


def _relative_path(path: str, *, output_dir: Path, repo_root: Path) -> str | None:
    if not path:
        return None
    candidate = Path(path)
    try:
        resolved = candidate.resolve()
    except Exception:
        resolved = candidate
    for base in (output_dir, repo_root):
        try:
            return str(resolved.relative_to(base.resolve()))
        except Exception:
            pass
    return candidate.name if candidate.name else None


def _reproducibility_role(role: str, generated_by: str | None) -> str:
    if generated_by:
        if role in {"selected_program", "metadata_json", "gmt", "run_summary", "full_program"}:
            return "final_output"
        return "intermediate_output"
    if role.startswith("resource:") or role in {
        "gtf",
        "sample_metadata_tsv",
        "comparisons_tsv",
        "subject_metadata_tsv",
        "cell_metadata_tsv",
        "counts_tsv",
        "deg_tsv",
        "feature_mapping_tsv",
    }:
        return "required_input"
    return "supporting_input"


def _checksum_records(path: Path | None, node: dict[str, Any]) -> list[dict[str, str]]:
    c2m2 = node.get("c2m2_properties", {})
    md5 = c2m2.get("md5") if isinstance(c2m2, dict) else None
    if path is not None and path.exists() and path.is_file():
        return checksum_list(path, str(md5) if isinstance(md5, str) and md5 else None)
    if isinstance(md5, str) and md5:
        return [{"algorithm": "md5", "encoding": "base64", "value": md5}]
    return []


def _extract_paths(metadata_payload: dict[str, Any], compact_graph: dict[str, Any]) -> dict[tuple[str, str], str]:
    paths: dict[tuple[str, str], str] = {}
    for section in ("input", "output"):
        content = metadata_payload.get(section, {})
        if not isinstance(content, dict):
            continue
        for item in content.get("files", []):
            if not isinstance(item, dict):
                continue
            role = str(item.get("role", "")).strip()
            path = str(item.get("local_path") or item.get("path") or "").strip()
            if role and path:
                paths[(role, Path(path).name)] = path
    for node in compact_graph.get("nodes", []):
        if not isinstance(node, dict) or node.get("type") != "File":
            continue
        role = _role_from_compact_node(node)
        c2m2 = node.get("c2m2_properties", {})
        if not isinstance(c2m2, dict):
            continue
        local_id = str(c2m2.get("local_id", "")).strip()
        filename = str(c2m2.get("filename", "")).strip()
        if role and local_id:
            paths[(role, filename or Path(local_id).name)] = local_id
    return paths


def _step_order(activities: dict[str, dict[str, Any]], edges: list[dict[str, Any]]) -> list[str]:
    deps: dict[str, set[str]] = {node_id: set() for node_id in activities}
    producer_by_file: dict[str, str] = {}
    consumers_by_file: dict[str, set[str]] = {}
    for edge in edges:
        source = str(edge.get("source", ""))
        target = str(edge.get("target", ""))
        if source in activities:
            producer_by_file[target] = source
        if target in activities:
            consumers_by_file.setdefault(source, set()).add(target)
    for file_id, producer in producer_by_file.items():
        for consumer in consumers_by_file.get(file_id, set()):
            deps.setdefault(consumer, set()).add(producer)
    ordered: list[str] = []
    remaining = {key: set(value) for key, value in deps.items()}
    while remaining:
        ready = sorted([key for key, value in remaining.items() if not value])
        if not ready:
            ordered.extend(sorted(remaining))
            break
        for key in ready:
            ordered.append(key)
            remaining.pop(key, None)
        for value in remaining.values():
            value.difference_update(ready)
    return ordered


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


def _entrypoint_to_base_argv(entrypoint: str) -> list[str]:
    parts = entrypoint.strip().split()
    if len(parts) >= 3 and parts[0] == "geneset-extractors":
        return [sys.executable, "-m", "geneset_extractors.cli", *parts[1:]]
    return [sys.executable, "-m", "geneset_extractors.cli"]


def _reconstruct_step_command(entrypoint: str, parameters: dict[str, Any], used_files: list[dict[str, Any]]) -> str:
    argv = _entrypoint_to_base_argv(entrypoint)
    for item in sorted(used_files, key=lambda x: str(x.get("role", ""))):
        role = str(item.get("role", "")).strip()
        path = str(item.get("path", "")).strip()
        if role and path:
            argv.extend([f"--{role}", path])
    for key in sorted(parameters):
        _append_parameter_tokens(argv, key, parameters[key])
    return shlex.join(argv)


def _portable_command(command: str, file_records: list[dict[str, Any]]) -> str:
    portable = command
    replacements: list[tuple[str, str]] = []
    for item in file_records:
        absolute = str(item.get("absolute_path_observed") or "").strip()
        relative = str(item.get("path") or "").strip()
        if absolute and relative:
            replacements.append((absolute, relative))
    for absolute, relative in sorted(replacements, key=lambda item: len(item[0]), reverse=True):
        portable = portable.replace(absolute, relative)
    return portable


def _command_matches_entrypoint(command: str, entrypoint: str) -> bool:
    required = [token for token in entrypoint.split() if token != "geneset-extractors"]
    lowered = command.lower()
    return all(token.lower() in lowered for token in required)


def _classify_input_edge(role: str) -> str:
    lowered = role.lower()
    metadata_tokens = ("gtf", "alias", "annotation", "metadata", "reference", "resource", "template", "manifest", "ontology")
    if any(token in lowered for token in metadata_tokens):
        return "metadata input"
    return "data input"


def _build_replay_file(
    *,
    node: dict[str, Any],
    generated_by: str | None,
    output_dir: Path,
    repo_root: Path,
    overlay: dict[str, Any],
    known_paths: dict[tuple[str, str], str],
) -> dict[str, Any]:
    role = _role_from_compact_node(node)
    c2m2 = node.get("c2m2_properties", {}) if isinstance(node.get("c2m2_properties"), dict) else {}
    filename = str(c2m2.get("filename") or node.get("name") or role)
    local_path = str(c2m2.get("local_id") or known_paths.get((role, filename), "")).strip()
    path_obj = Path(local_path) if local_path else None
    relative = _relative_path(local_path, output_dir=output_dir, repo_root=repo_root) if local_path else None
    source = source_for_role(overlay, role)
    mirror = mirror_for_role(overlay, role)
    return {
        "id": f"#{node['id']}",
        "name": str(node.get("name", filename)),
        "role": role,
        "reproducibility_role": _reproducibility_role(role, generated_by),
        "path": relative,
        "absolute_path_observed": local_path or None,
        "size_in_bytes": file_size(path_obj) if path_obj is not None and path_obj.exists() and path_obj.is_file() else int(c2m2.get("size_in_bytes", 0) or 0),
        "checksums": _checksum_records(path_obj, node),
        "description": str(node.get("description", "")),
        "source": source,
        "published": mirror,
        "c2m2": c2m2 if c2m2 else None,
        "dcc_url": node.get("dcc_url"),
        "drc_url": node.get("drc_url"),
        "generated_by": generated_by,
    }


def _materialize_environment_outputs(output_dir: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    files: list[dict[str, Any]] = []
    software_env: dict[str, Any] = {}
    pip_text = pip_freeze_text()
    if pip_text:
        pip_path = output_dir / "provenance.pip.freeze.txt"
        pip_path.write_text(pip_text + ("\n" if not pip_text.endswith("\n") else ""), encoding="utf-8")
        files.append(
            {
                "id": "#file:provenance_pip_freeze",
                "name": pip_path.name,
                "role": "provenance_pip_freeze",
                "reproducibility_role": "supporting_output",
                "path": pip_path.name,
                "absolute_path_observed": str(pip_path),
                "size_in_bytes": file_size(pip_path),
                "checksums": checksum_list(pip_path),
                "description": "Captured pip freeze for replay environment audit.",
                "source": None,
                "published": None,
                "c2m2": None,
                "dcc_url": pip_path.resolve().as_uri(),
                "drc_url": pip_path.resolve().as_uri(),
                "generated_by": None,
            }
        )
        software_env["pip_freeze_path"] = pip_path.name
    r_text = r_session_info_text()
    if r_text:
        r_path = output_dir / "provenance.r_session_info.txt"
        r_path.write_text(r_text + ("\n" if not r_text.endswith("\n") else ""), encoding="utf-8")
        files.append(
            {
                "id": "#file:provenance_r_session_info",
                "name": r_path.name,
                "role": "provenance_r_session_info",
                "reproducibility_role": "supporting_output",
                "path": r_path.name,
                "absolute_path_observed": str(r_path),
                "size_in_bytes": file_size(r_path),
                "checksums": checksum_list(r_path),
                "description": "Captured R sessionInfo for replay environment audit.",
                "source": None,
                "published": None,
                "c2m2": None,
                "dcc_url": r_path.resolve().as_uri(),
                "drc_url": r_path.resolve().as_uri(),
                "generated_by": None,
            }
        )
        software_env["r_session_info_path"] = r_path.name
    return files, software_env


def _derive_linked_data(replay: dict[str, Any]) -> dict[str, Any]:
    graph: list[dict[str, Any]] = []
    graph.append(
        {
            "@id": replay["geneset"]["id"],
            "@type": ["prov:Entity", "dig:GeneSet"],
            "name": replay["geneset"]["name"],
            "description": replay["geneset"].get("description"),
        }
    )
    software = replay.get("software", {})
    repository = software.get("repository")
    if isinstance(repository, dict):
        graph.append(
            {
                "@id": repository["id"],
                "@type": ["prov:Entity", "dig:SoftwareRepository"],
                "dig:repo_url": repository.get("repo_url"),
                "dig:git_commit": repository.get("git_commit"),
                "dig:git_ref": repository.get("git_ref"),
                "dig:git_dirty": repository.get("git_dirty"),
            }
        )
    environment = software.get("environment")
    if isinstance(environment, dict):
        graph.append(
            {
                "@id": environment["id"],
                "@type": ["prov:Entity", "dig:SoftwareEnvironment"],
                "platform": environment.get("platform"),
                "dig:python_version": environment.get("python_version"),
                "dig:r_version": environment.get("r_version"),
            }
        )
    all_files = list(replay.get("inputs", []))
    for step in replay.get("steps", []):
        for item in step.get("generated_files", []):
            all_files.append(item)
    seen_file_ids: set[str] = set()
    for item in all_files:
        if not isinstance(item, dict):
            continue
        file_id = item["id"]
        if file_id in seen_file_ids:
            continue
        seen_file_ids.add(file_id)
        graph.append(
            {
                "@id": file_id,
                "@type": ["prov:Entity", "schema:MediaObject", "dig:LocalAnalysisFile"],
                "name": item.get("name"),
                "description": item.get("description"),
                "dig:file_role": item.get("role"),
                "dig:reproducibility_role": item.get("reproducibility_role"),
                "dig:relative_path": item.get("path"),
                "dig:absolute_path_observed": item.get("absolute_path_observed"),
                "contentSize": item.get("size_in_bytes"),
                "dig:checksums": item.get("checksums", []),
            }
        )
        source = item.get("source")
        if isinstance(source, dict) and source.get("remote_url"):
            graph.append(
                {
                    "@id": f"#remote:{item['role']}",
                    "@type": ["prov:Entity", "schema:DataDownload", "dig:RemoteSourceFile"],
                    "name": source.get("name", f"{item['role']} source"),
                    "contentUrl": source.get("remote_url"),
                    "dig:download_command": source.get("download_command"),
                }
            )
        published = item.get("published")
        if isinstance(published, dict) and (published.get("remote_url") or published.get("s3_url") or published.get("uri")):
            graph.append(
                {
                    "@id": f"#published:{item['role']}",
                    "@type": ["prov:Entity", "schema:DataDownload", "dig:PublishedMirror"],
                    "name": published.get("name", f"{item['role']} published"),
                    "contentUrl": published.get("remote_url") or published.get("s3_url") or published.get("uri"),
                }
            )
    for step in replay.get("steps", []):
        graph.append(
            {
                "@id": step["id"],
                "@type": ["prov:Activity", "dig:CommandStep"],
                "name": step.get("name"),
                "dig:step_index": step.get("step_index"),
                "dig:entrypoint": step.get("entrypoint"),
                "dig:observed_command": step.get("observed_command"),
                "dig:replay_command": step.get("replay_command"),
                "dig:command_argv": step.get("command_argv"),
                "prov:used": [{"@id": item} for item in step.get("used", [])],
                "prov:generated": [{"@id": item["id"]} for item in step.get("generated_files", [])],
                "prov:wasAssociatedWith": [{"@id": repository["id"]}, {"@id": environment["id"]}] if isinstance(repository, dict) and isinstance(environment, dict) else [],
            }
        )
    return {"@context": CONTEXT, "@graph": graph}


def _derive_c2m2(replay: dict[str, Any]) -> dict[str, Any]:
    nodes: list[dict[str, Any]] = []
    edges: list[dict[str, Any]] = []
    geneset = replay["geneset"]
    geneset_id = str(geneset["id"]).lstrip("#")
    nodes.append(
        {
            "id": geneset_id,
            "type": "GeneSet",
            "name": geneset.get("name"),
            "description": geneset.get("description"),
            "dcc_url": geneset.get("dcc_url") or f"urn:geneset:{geneset_id}",
            "drc_url": geneset.get("drc_url") or geneset.get("dcc_url") or f"urn:geneset:{geneset_id}",
            "c2m2_properties": {
                "id": geneset_id,
                "name": geneset.get("name"),
                "description": geneset.get("description"),
            },
        }
    )
    all_files = list(replay.get("inputs", []))
    for step in replay.get("steps", []):
        for item in step.get("generated_files", []):
            all_files.append(item)
    seen_files: set[str] = set()
    for item in all_files:
        if not isinstance(item, dict):
            continue
        file_id = str(item["id"]).lstrip("#")
        if file_id in seen_files:
            continue
        seen_files.add(file_id)
        c2m2 = item.get("c2m2") or {}
        nodes.append(
            {
                "id": file_id,
                "type": "File",
                "name": item.get("name"),
                "description": item.get("description"),
                "dcc_url": item.get("dcc_url") or item.get("absolute_path_observed") or f"urn:file:{file_id}",
                "drc_url": item.get("drc_url") or item.get("absolute_path_observed") or f"urn:file:{file_id}",
                "c2m2_properties": {
                    "filename": c2m2.get("filename") or item.get("name"),
                    "persistent_id": c2m2.get("persistent_id") or file_id,
                    "local_id": c2m2.get("local_id") or item.get("absolute_path_observed") or item.get("path") or item.get("name"),
                    "size_in_bytes": c2m2.get("size_in_bytes") or item.get("size_in_bytes") or 0,
                    "_uuid": c2m2.get("_uuid") or file_id,
                    "md5": next((c["value"] for c in item.get("checksums", []) if c.get("algorithm") == "md5"), ""),
                },
            }
        )
    for step in replay.get("steps", []):
        step_id = str(step["id"]).lstrip("#")
        nodes.append(
            {
                "id": step_id,
                "type": "AnalysisType",
                "name": step.get("name"),
                "description": step.get("description") or f"Replay step {step.get('name')}",
                "dcc_url": replay["software"]["repository"].get("repo_url", ""),
                "drc_url": replay["software"]["repository"].get("repo_url", ""),
                "c2m2_properties": {
                    "id": step_id,
                    "name": step.get("name"),
                    "description": step.get("description") or f"Replay step {step.get('name')}",
                    "synonyms": [],
                },
                "analysis": {
                    "script_url": step.get("script_url") or replay["software"]["repository"].get("repo_url", ""),
                    "command": step.get("observed_command"),
                    "parameters": step.get("parameters", {}),
                    "environment": {
                        "entrypoint": step.get("entrypoint"),
                        "container_image": step.get("container_image"),
                    },
                },
            }
        )
        for used_id in step.get("used", []):
            role = ""
            for item in replay.get("inputs", []):
                if item.get("id") == used_id:
                    role = str(item.get("role", ""))
                    break
            edges.append(
                {
                    "id": f"{used_id.lstrip('#')}_to_{step_id}",
                    "source": used_id.lstrip("#"),
                    "target": step_id,
                    "label": _classify_input_edge(role),
                    "description": f"{role or 'input'} used by {step.get('name')}",
                }
            )
        for item in step.get("generated_files", []):
            edges.append(
                {
                    "id": f"{step_id}_to_{str(item['id']).lstrip('#')}",
                    "source": step_id,
                    "target": str(item["id"]).lstrip("#"),
                    "label": "data output",
                    "description": f"{step.get('name')} generated {item.get('role')}",
                }
            )
        edges.append(
            {
                "id": f"{step_id}_to_{geneset_id}",
                "source": step_id,
                "target": geneset_id,
                "label": "data output",
                "description": f"{step.get('name')} contributed to gene set output",
            }
        )
    return {"nodes": nodes, "edges": edges}


def build_single_record(
    *,
    metadata_payload: dict[str, Any],
    compact_payload: dict[str, Any],
    overlay: dict[str, Any],
    output_dir: Path,
    repo_root: Path,
) -> dict[str, Any]:
    graph_key = next(iter(compact_payload))
    compact_graph = compact_payload[graph_key]
    nodes = [node for node in compact_graph.get("nodes", []) if isinstance(node, dict)]
    edges = [edge for edge in compact_graph.get("edges", []) if isinstance(edge, dict)]
    known_paths = _extract_paths(metadata_payload, compact_graph)

    activities = {str(node["id"]): node for node in nodes if node.get("type") == "AnalysisType"}
    generated_by_file: dict[str, str] = {}
    used_by_step: dict[str, list[str]] = {activity_id: [] for activity_id in activities}
    generated_by_step: dict[str, list[str]] = {activity_id: [] for activity_id in activities}
    file_by_id: dict[str, dict[str, Any]] = {}
    for node in nodes:
        if node.get("type") == "File":
            file_by_id[str(node["id"])] = node
    for edge in edges:
        source = str(edge.get("source", ""))
        target = str(edge.get("target", ""))
        if source in activities:
            generated_by_file[target] = source
            generated_by_step[source].append(target)
        if target in activities:
            used_by_step[target].append(source)

    replay_inputs: list[dict[str, Any]] = []
    generated_files_by_step: dict[str, list[dict[str, Any]]] = {activity_id: [] for activity_id in activities}
    all_replay_files: dict[str, dict[str, Any]] = {}
    for file_id, node in file_by_id.items():
        generated_by = generated_by_file.get(file_id)
        replay_file = _build_replay_file(
            node=node,
            generated_by=f"#{generated_by}" if generated_by else None,
            output_dir=output_dir,
            repo_root=repo_root,
            overlay=overlay,
            known_paths=known_paths,
        )
        all_replay_files[replay_file["id"]] = replay_file
        if generated_by:
            generated_files_by_step[generated_by].append(replay_file)
        else:
            replay_inputs.append(replay_file)

    support_files, env_file_bits = _materialize_environment_outputs(output_dir)
    for item in support_files:
        all_replay_files[item["id"]] = item

    order = _step_order(activities, edges)
    steps: list[dict[str, Any]] = []
    for index, activity_id in enumerate(order, start=1):
        node = activities[activity_id]
        analysis = node.get("analysis", {}) if isinstance(node.get("analysis"), dict) else {}
        env = analysis.get("environment", {}) if isinstance(analysis.get("environment"), dict) else {}
        entrypoint = str(env.get("entrypoint") or node.get("name") or activity_id)
        parameters = analysis.get("parameters") if isinstance(analysis.get("parameters"), dict) else {}
        used_files = [all_replay_files[f"#{item}"] for item in used_by_step.get(activity_id, []) if f"#{item}" in all_replay_files]
        reconstructed = _reconstruct_step_command(entrypoint, parameters, used_files)
        raw_command = str(analysis.get("command") or "")
        observed_command = raw_command if raw_command and _command_matches_entrypoint(raw_command, entrypoint) else reconstructed
        replay_command = _portable_command(reconstructed, list(all_replay_files.values()))
        steps.append(
            {
                "id": f"#{activity_id}",
                "name": str(node.get("name", activity_id)),
                "description": str(node.get("description", "")),
                "step_index": index,
                "entrypoint": entrypoint,
                "observed_command": observed_command,
                "replay_command": replay_command,
                "command_argv": shlex.split(replay_command),
                "parameters": parameters,
                "parameter_fingerprint": stable_hash_object(parameters),
                "script_url": analysis.get("script_url"),
                "container_image": env.get("container_image"),
                "used": [item["id"] for item in used_files],
                "generated_files": generated_files_by_step.get(activity_id, []),
            }
        )

    repo_info = git_info(repo_root)
    software = software_overlay(overlay) or {}
    py_env = python_environment()
    software_section = {
        "repository": {
            "id": "#repo:dig_gene_set_extractors",
            "repo_url": str(software.get("repo_url") or metadata_payload.get("converter", {}).get("code", {}).get("repo_url") or ""),
            "git_commit": str(software.get("git_commit") or repo_info.get("git_commit") or "unknown"),
            "git_ref": software.get("git_ref") or repo_info.get("git_ref"),
            "git_dirty": bool(software.get("git_dirty")) if "git_dirty" in software else bool(repo_info.get("git_dirty")),
            "install_command": software.get("install_command"),
        },
        "environment": {
            "id": "#env:rnaseq",
            "platform": str(py_env.get("platform")),
            "python_version": str(py_env.get("python_version")),
            "r_version": r_package_version("base"),
            **env_file_bits,
        },
    }

    geneset_id = f"#geneset:{metadata_payload.get('geneset_id', graph_key)}"
    geneset = {
        "id": geneset_id,
        "focus_node_id": metadata_payload.get("provenance", {}).get("focus_node_id"),
        "geneset_id": metadata_payload.get("geneset_id", graph_key),
        "name": metadata_payload.get("gene_set", {}).get("name", graph_key),
        "description": metadata_payload.get("gene_set", {}).get("description", ""),
        "signature_name": metadata_payload.get("converter", {}).get("parameters", {}).get("signature_name"),
        "comparison_id": metadata_payload.get("converter", {}).get("parameters", {}).get("comparison_label"),
        "tissue_id": metadata_payload.get("converter", {}).get("parameters", {}).get("dataset_label"),
        "dataset": metadata_payload.get("converter", {}).get("parameters", {}).get("dataset_label"),
        "dataset_version": metadata_payload.get("standard_version"),
        "model_fingerprint": stable_hash_object(metadata_payload.get("converter", {}).get("parameters", {})),
    }

    target_output = None
    final_files = [item for item in all_replay_files.values() if item.get("reproducibility_role") == "final_output"]
    for preferred_role in ("gmt", "selected_program"):
        target_output = next((item for item in final_files if item.get("role") == preferred_role), None)
        if target_output:
            break
    if target_output is None and final_files:
        target_output = final_files[0]

    replay = {
        "geneset": geneset,
        "software": software_section,
        "inputs": replay_inputs + support_files,
        "steps": steps,
        "target_output": {
            "id": target_output["id"] if target_output else None,
            "role": target_output.get("role") if target_output else None,
            "path": target_output.get("path") if target_output else None,
            "sha256": next((c["value"] for c in target_output.get("checksums", []) if c.get("algorithm") == "sha256"), None) if target_output else None,
            "published": target_output.get("published") if target_output else None,
        },
    }

    return {
        "schema_version": "geneset-provenance/1.0",
        "record_type": "single_geneset",
        "replay": replay,
        "linked_data": _derive_linked_data(replay),
        "c2m2": _derive_c2m2(replay),
    }
