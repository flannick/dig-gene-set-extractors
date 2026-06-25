from __future__ import annotations

import json
from pathlib import Path
from string import Formatter
from typing import Any

from geneset_extractors.core.metadata import write_provenance_from_metadata
from geneset_extractors.core.provenance import write_canonical_json


def load_metadata(path: str | Path) -> dict[str, Any]:
    meta_path = Path(path)
    payload = json.loads(meta_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("metadata payload must be a JSON object")
    return payload


def load_model_sidecar(metadata_path: str | Path) -> dict[str, Any] | None:
    meta_path = Path(metadata_path)
    sidecar_path = meta_path.with_name("geneset.model.json")
    if not sidecar_path.exists():
        return None
    payload = json.loads(sidecar_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("model sidecar payload must be a JSON object")
    return payload


def flatten_template_context(payload: dict[str, Any], model_payload: dict[str, Any] | None = None) -> dict[str, str]:
    context: dict[str, str] = {}

    def add_flat(prefix: str, value: Any) -> None:
        if isinstance(value, dict):
            for key, inner in value.items():
                child_prefix = f"{prefix}.{key}" if prefix else str(key)
                add_flat(child_prefix, inner)
            return
        if isinstance(value, list):
            return
        if not prefix:
            return
        context[prefix] = "" if value is None else str(value)

    converter = payload.get("converter", {})
    if isinstance(converter, dict):
        parameters = converter.get("parameters", {})
        if isinstance(parameters, dict):
            for key, value in parameters.items():
                if not isinstance(value, (dict, list)):
                    context[str(key)] = "" if value is None else str(value)

    gene_set = payload.get("gene_set", {})
    if isinstance(gene_set, dict):
        for key in ("id", "name", "description", "assay", "data_type", "organism", "genome_build"):
            value = gene_set.get(key)
            context[f"gene_set_{key}"] = "" if value is None else str(value)

    for key in ("geneset_id", "standard_name", "standard_version", "schema_version"):
        value = payload.get(key)
        context[str(key)] = "" if value is None else str(value)

    add_flat("", payload)
    if isinstance(model_payload, dict):
        add_flat("model", model_payload)
    return dict(sorted(context.items()))


def build_template_context(metadata_path: str | Path, payload: dict[str, Any] | None = None) -> dict[str, str]:
    meta_path = Path(metadata_path)
    metadata_payload = payload if payload is not None else load_metadata(meta_path)
    model_payload = load_model_sidecar(meta_path)
    return flatten_template_context(metadata_payload, model_payload)


def render_template(template: str, context: dict[str, str]) -> str:
    missing = []
    for _, field_name, _, _ in Formatter().parse(template):
        if not field_name:
            continue
        if field_name not in context:
            missing.append(field_name)
    if missing:
        raise ValueError("Template references unknown variable(s): {0}".format(", ".join(sorted(set(missing)))))
    rendered_parts: list[str] = []
    for literal_text, field_name, format_spec, conversion in Formatter().parse(template):
        if literal_text:
            rendered_parts.append(literal_text)
        if field_name is None:
            continue
        value = context[field_name]
        if conversion:
            if conversion == "r":
                value = repr(value)
            elif conversion == "s":
                value = str(value)
            elif conversion == "a":
                value = ascii(value)
            else:
                raise ValueError(f"Unsupported template conversion: !{conversion}")
        if format_spec:
            value = format(value, format_spec)
        rendered_parts.append(str(value))
    return "".join(rendered_parts)


def set_dotted_value(payload: dict[str, Any], dotted_key: str, value: str) -> None:
    parts = [part for part in str(dotted_key).split(".") if part]
    if not parts:
        raise ValueError("Empty dotted key is not allowed")
    current: dict[str, Any] = payload
    for part in parts[:-1]:
        existing = current.get(part)
        if existing is None:
            current[part] = {}
            existing = current[part]
        if not isinstance(existing, dict):
            raise ValueError("Cannot descend into non-object field: {0}".format(".".join(parts[:-1])))
        current = existing
    current[parts[-1]] = value


def apply_metadata_patch(
    *,
    metadata_path: str | Path,
    meta_out: str | Path | None = None,
    provenance_out: str | Path | None = None,
    description_template: str | None = None,
    gene_set_description: str | None = None,
    set_values: list[tuple[str, str]] | None = None,
    provenance_overlay_json: str | None = None,
    upstream_provenance_graph_path: str | None = None,
    provenance_mirror_local_prefix: str | None = None,
    provenance_mirror_remote_prefix: str | None = None,
) -> dict[str, str]:
    meta_path = Path(metadata_path)
    payload = load_metadata(meta_path)

    if description_template:
        context = build_template_context(meta_path, payload)
        set_dotted_value(payload, "gene_set.description", render_template(description_template, context))
    if gene_set_description is not None:
        set_dotted_value(payload, "gene_set.description", gene_set_description)
    for dotted_key, value in (set_values or []):
        set_dotted_value(payload, dotted_key, value)

    out_meta_path = Path(meta_out) if meta_out is not None else meta_path
    write_canonical_json(out_meta_path, payload)
    out_prov_path = write_provenance_from_metadata(
        out_meta_path,
        provenance_path=provenance_out,
        provenance_overlay_json=provenance_overlay_json,
        upstream_provenance_graph_path=upstream_provenance_graph_path,
        provenance_mirror_local_prefix=provenance_mirror_local_prefix,
        provenance_mirror_remote_prefix=provenance_mirror_remote_prefix,
    )
    return {
        "metadata_path": str(out_meta_path),
        "provenance_path": str(out_prov_path),
    }
