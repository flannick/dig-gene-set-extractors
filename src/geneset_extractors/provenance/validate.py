from __future__ import annotations

from typing import Any


def validate_single_record(payload: dict[str, Any]) -> None:
    required = {"schema_version", "record_type", "replay", "linked_data", "c2m2"}
    missing = sorted(required - set(payload))
    if missing:
        raise ValueError("missing required provenance keys: " + ", ".join(missing))
    if payload.get("record_type") != "single_geneset":
        raise ValueError("record_type must be 'single_geneset'")
    replay = payload.get("replay")
    if not isinstance(replay, dict):
        raise ValueError("replay must be an object")
    for key in ("geneset", "software", "inputs", "steps", "target_output"):
        if key not in replay:
            raise ValueError(f"replay missing required key {key!r}")
    if not isinstance(replay["inputs"], list):
        raise ValueError("replay.inputs must be a list")
    if not isinstance(replay["steps"], list) or not replay["steps"]:
        raise ValueError("replay.steps must be a non-empty list")
    linked = payload.get("linked_data")
    if not isinstance(linked, dict) or not isinstance(linked.get("@graph"), list):
        raise ValueError("linked_data must contain @graph array")
    c2m2 = payload.get("c2m2")
    if not isinstance(c2m2, dict):
        raise ValueError("c2m2 must be an object")
    if not isinstance(c2m2.get("nodes"), list) or not isinstance(c2m2.get("edges"), list):
        raise ValueError("c2m2 must contain nodes and edges arrays")
