from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .model import GeneSetProvenanceRecord


CONTEXT = {
    "schema": "https://schema.org/",
    "prov": "http://www.w3.org/ns/prov#",
    "dig": "https://dig.example.org/provenance#",
    "name": "schema:name",
    "description": "schema:description",
    "contentUrl": "schema:contentUrl",
    "contentSize": "schema:contentSize",
    "sha256": "schema:sha256",
}


def node_to_jsonld(node: Any) -> dict[str, Any]:
    if isinstance(node, dict):
        return node
    if hasattr(node, "to_jsonld"):
        return node.to_jsonld()
    raise TypeError(f"Cannot serialize node of type {type(node)}")


def to_jsonld(record: GeneSetProvenanceRecord) -> dict[str, Any]:
    return {
        "@context": CONTEXT,
        "@id": record.record_id,
        "@type": ["dig:GeneSetProvenanceRecord", "prov:Entity"],
        "schema_version": record.schema_version,
        "geneset": record.geneset,
        "model": record.model,
        "replay_plan": {"@id": record.replay_plan_id},
        "@graph": [node_to_jsonld(node) for node in record.graph],
    }


def validate_references(payload: dict[str, Any]) -> None:
    ids = {payload.get("@id")}
    ids |= {node.get("@id") for node in payload.get("@graph", []) if isinstance(node, dict)}
    missing: list[str] = []

    def walk(obj: Any) -> None:
        if isinstance(obj, dict):
            if set(obj.keys()) == {"@id"} and obj["@id"] not in ids:
                missing.append(str(obj["@id"]))
            for value in obj.values():
                walk(value)
        elif isinstance(obj, list):
            for value in obj:
                walk(value)

    walk(payload)
    if missing:
        raise ValueError("Unresolved provenance references: " + ", ".join(sorted(set(missing))))


def write_jsonld(record: GeneSetProvenanceRecord, path: Path) -> None:
    payload = to_jsonld(record)
    validate_references(payload)
    path.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")

