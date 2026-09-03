"""Write DAPPER-ID-1 provenance YAML from DIG's legacy provenance graph.

The structural conversion is adapted from DAPPER's
``schema/converter/geneset_to_dapper.py``. The identifier implementation is
the DAPPER-ID-1 algorithm for the three DAPPER node classes emitted here,
pinned to DAPPER schema revision ``57331c2877b223a8dd691266c877e173049d65f3``.
"""
from __future__ import annotations

import base64
import hashlib
import json
from pathlib import Path
from typing import Any

import yaml
from rdflib import Graph, Literal, URIRef
from rdflib.compare import to_canonical_graph


_INPUT_LABELS = {"data input", "metadata input"}
_EDGE_ROLE = {"data input": "data_input", "metadata input": "metadata_input"}
_OUTPUT_LABEL = "data output"
_BUCKETS = (
    "c2m2_files",
    "activities",
    "gene_sets",
    "used_edges",
    "was_generated_by_edges",
)
DAPPER_SCHEMA_REVISION = "57331c2877b223a8dd691266c877e173049d65f3"
_CLASS_BY_BUCKET = {"c2m2_files": "C2M2File", "activities": "Activity", "gene_sets": "GeneSet"}
_HASHABLE_SLOTS = {
    "C2M2File": ("c2m2_uuid", "dcc_url", "description", "drc_url", "filename", "local_id", "md5", "name", "persistent_id", "sha256", "size_in_bytes"),
    "Activity": ("activity_type", "code_version", "command", "container_image", "dcc_url", "description", "drc_url", "entrypoint", "has_agentic_workspace", "has_lineage_step", "name", "observed_command", "repo_url", "script_url", "software_name", "software_version"),
    "GeneSet": ("access_level", "alternate_identifier", "assay", "controlled_access", "data_type", "dcc_url", "description", "drc_url", "funded_by", "genome_build", "has_contributor", "has_creator", "has_license", "has_recommended_citation", "is_described_by", "member_type", "members", "n_genes", "n_members", "n_sets", "name", "organism", "term", "term_prefix", "was_attributed_to", "was_derived_from", "was_generated_by"),
}
_SELF = URIRef("urn:dapper:self")
_CLASS_PREDICATE = URIRef("urn:dapper:class")
_SLOT_NAMESPACE = "urn:dapper:slot:"
_SELF_REFERENCE = "urn:dapper:self-reference"


def _sha512t24u(payload: bytes) -> str:
    return base64.urlsafe_b64encode(hashlib.sha512(payload).digest()[:24]).decode("ascii")


def _dapper_digest(value: Any) -> str | None:
    if not isinstance(value, str) or not value.startswith("dapper:"):
        return None
    _, _, suffix = value.partition(":")
    _, dot, digest = suffix.partition(".")
    return digest if dot and digest else None


def _replace_self(value: Any, identifier: str | None) -> Any:
    if not identifier:
        return value
    if isinstance(value, str):
        return value.replace(identifier, _SELF_REFERENCE)
    if isinstance(value, list):
        return [_replace_self(item, identifier) for item in value]
    if isinstance(value, dict):
        return {key: _replace_self(item, identifier) for key, item in value.items()}
    return value


def _literal(value: Any) -> Literal:
    digest = _dapper_digest(value)
    if digest is not None:
        value = digest
    if isinstance(value, (dict, list)):
        value = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return Literal(value)


def _compute_id(node: dict[str, Any], class_name: str, identifier: str | None) -> str:
    """Compute DAPPER-ID-1 exactly for DIG's emitted DAPPER node classes."""
    graph = Graph()
    graph.add((_SELF, _CLASS_PREDICATE, Literal(class_name)))
    for slot in _HASHABLE_SLOTS[class_name]:
        value = _replace_self(node.get(slot), identifier)
        if value is None:
            continue
        if isinstance(value, list):
            for index, item in enumerate(value):
                if item is not None:
                    graph.add((_SELF, URIRef(f"{_SLOT_NAMESPACE}{slot}[{index}]"), _literal(item)))
        else:
            graph.add((_SELF, URIRef(f"{_SLOT_NAMESPACE}{slot}"), _literal(value)))
    canonical = to_canonical_graph(graph)
    lines = sorted(line for line in canonical.serialize(format="nt").split("\n") if line.strip())
    digest = _sha512t24u("\n".join(lines).encode("utf-8"))
    return f"dapper:{class_name}.{digest}"


def _rewrite(value: Any, identifiers: dict[str, str]) -> Any:
    if isinstance(value, str):
        for old, new in sorted(identifiers.items(), key=lambda item: len(item[0]), reverse=True):
            value = value.replace(old, new)
        return value
    if isinstance(value, list):
        return [_rewrite(item, identifiers) for item in value]
    if isinstance(value, dict):
        return {key: _rewrite(item, identifiers) for key, item in value.items()}
    return value


def _rewrite_node(node: dict[str, Any], identifiers: dict[str, str]) -> dict[str, Any]:
    """Mirror DAPPER's rewrite rule for the literal-only emitted node fields."""
    rewritten: dict[str, Any] = {}
    for key, value in node.items():
        if key == "id":
            rewritten[key] = identifiers.get(value, value)
        elif isinstance(value, str) and (" " in value or "\n" in value):
            rewritten[key] = _rewrite(value, identifiers)
        elif isinstance(value, (list, dict)):
            # DIG's emitted DAPPER node slots are literals, not relationships;
            # list/dict members therefore use exact rather than substring rewrites.
            rewritten[key] = _rewrite_exact(value, identifiers)
        else:
            rewritten[key] = value
    return rewritten


def _rewrite_exact(value: Any, identifiers: dict[str, str]) -> Any:
    if isinstance(value, str):
        return identifiers.get(value, value)
    if isinstance(value, list):
        return [_rewrite_exact(item, identifiers) for item in value]
    if isinstance(value, dict):
        return {key: _rewrite_exact(item, identifiers) for key, item in value.items()}
    return value


def _mint_dapper_ids(document: dict[str, list[dict[str, Any]]]) -> None:
    """Replace DIG node IDs and every edge reference with DAPPER-ID-1 IDs."""
    identifiers: dict[str, str] = {}
    for bucket, class_name in _CLASS_BY_BUCKET.items():
        for node in document.get(bucket, []):
            old_id = node.get("id")
            if isinstance(old_id, str):
                identifiers[old_id] = _compute_id(node, class_name, old_id)
    for bucket in _CLASS_BY_BUCKET:
        document[bucket] = [_rewrite_node(node, identifiers) for node in document.get(bucket, [])]
    for bucket in ("used_edges", "was_generated_by_edges"):
        document[bucket] = [_rewrite(edge, identifiers) for edge in document.get(bucket, [])]


def _clean(payload: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in payload.items() if value not in (None, "", [], {})}


def _sha256_by_local_id(metadata: dict[str, Any]) -> dict[str, str]:
    checksums: dict[str, str] = {}
    for record in (metadata.get("input", {}) or {}).get("files", []) or []:
        if not isinstance(record, dict) or not record.get("sha256"):
            continue
        for key in ("path", "local_path"):
            value = record.get(key)
            if value:
                checksums[str(value)] = str(record["sha256"])
    return checksums


def _file(node: dict[str, Any], checksums: dict[str, str]) -> dict[str, Any]:
    c2m2 = node.get("c2m2_properties", {}) or {}
    local_id = c2m2.get("local_id")
    return _clean(
        {
            "id": node.get("id"),
            "name": node.get("name"),
            "description": node.get("description"),
            "filename": c2m2.get("filename"),
            "persistent_id": c2m2.get("persistent_id"),
            "local_id": local_id,
            "c2m2_uuid": c2m2.get("_uuid"),
            "md5": c2m2.get("md5"),
            "sha256": checksums.get(str(local_id)),
            "size_in_bytes": c2m2.get("size_in_bytes"),
            "dcc_url": node.get("dcc_url"),
            "drc_url": node.get("drc_url"),
        }
    )


def _activity(node: dict[str, Any]) -> dict[str, Any]:
    analysis = node.get("analysis", {}) or {}
    environment = analysis.get("environment", {}) or {}
    c2m2 = node.get("c2m2_properties", {}) or {}
    result = _clean(
        {
            "id": node.get("id"),
            "name": node.get("name"),
            "description": node.get("description"),
            "command": analysis.get("command"),
            "observed_command": analysis.get("observed_command"),
            "script_url": analysis.get("script_url"),
            "repo_url": environment.get("repo_url"),
            "code_version": analysis.get("version"),
            "entrypoint": environment.get("entrypoint"),
            "container_image": environment.get("container_image"),
            "dcc_url": node.get("dcc_url"),
            "drc_url": node.get("drc_url"),
        }
    )
    if c2m2.get("synonyms"):
        result["aliases"] = list(c2m2["synonyms"])
    return result


def _gene_set(node: dict[str, Any], metadata: dict[str, Any]) -> dict[str, Any]:
    gene_set = metadata.get("gene_set", {}) or {}
    summary = metadata.get("summary", {}) or {}
    parameters = (metadata.get("converter", {}) or {}).get("parameters", {}) or {}
    return _clean(
        {
            "id": node.get("id"),
            "name": node.get("name"),
            "description": node.get("description") or gene_set.get("description"),
            "member_type": "gene",
            "assay": gene_set.get("assay"),
            "data_type": gene_set.get("data_type"),
            "organism": gene_set.get("organism"),
            "genome_build": gene_set.get("genome_build"),
            "n_genes": gene_set.get("n_genes") or summary.get("n_genes"),
            "n_sets": summary.get("n_sets_emitted"),
            "term_prefix": parameters.get("term_prefix"),
            "dcc_url": node.get("dcc_url"),
            "drc_url": node.get("drc_url"),
        }
    )


def _convert_graph(graph: dict[str, Any], metadata: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    """Map one DIG ``{nodes, edges}`` graph to DAPPER document buckets."""
    document: dict[str, list[dict[str, Any]]] = {bucket: [] for bucket in _BUCKETS}
    checksums = _sha256_by_local_id(metadata)
    for node in graph.get("nodes", []):
        if not isinstance(node, dict):
            continue
        match node.get("type"):
            case "File":
                document["c2m2_files"].append(_file(node, checksums))
            case "AnalysisType":
                document["activities"].append(_activity(node))
            case "GeneSet":
                document["gene_sets"].append(_gene_set(node, metadata))

    for edge in graph.get("edges", []):
        if not isinstance(edge, dict):
            continue
        label = edge.get("label")
        if label in _INPUT_LABELS:
            document["used_edges"].append(
                _clean(
                    {
                        "subject": edge.get("target"),
                        "predicate": "prov:used",
                        "object": edge.get("source"),
                        "edge_role": _EDGE_ROLE[str(label)],
                    }
                )
            )
        elif label == _OUTPUT_LABEL:
            document["was_generated_by_edges"].append(
                _clean(
                    {
                        "subject": edge.get("target"),
                        "predicate": "prov:wasGeneratedBy",
                        "object": edge.get("source"),
                    }
                )
            )
    _mint_dapper_ids(document)
    return {key: value for key, value in document.items() if value}


def build_dapper_provenance(
    legacy_payload: dict[str, Any], metadata: dict[str, Any]
) -> dict[str, list[dict[str, Any]]]:
    """Convert all graphs in a legacy DIG payload into one DAPPER YAML document."""
    document: dict[str, list[dict[str, Any]]] = {bucket: [] for bucket in _BUCKETS}
    for graph in legacy_payload.values():
        if not isinstance(graph, dict) or "nodes" not in graph:
            continue
        for bucket, values in _convert_graph(graph, metadata).items():
            document[bucket].extend(values)
    return {key: value for key, value in document.items() if value}


def write_dapper_provenance(
    path: str | Path, legacy_payload: dict[str, Any], metadata: dict[str, Any]
) -> Path:
    """Write a deterministic, readable DAPPER YAML sidecar."""
    output_path = Path(path)
    output_path.write_text(
        yaml.safe_dump(
            build_dapper_provenance(legacy_payload, metadata),
            sort_keys=False,
            default_flow_style=False,
            allow_unicode=True,
        ),
        encoding="utf-8",
    )
    return output_path
