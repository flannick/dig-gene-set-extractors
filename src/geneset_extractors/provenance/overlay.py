from __future__ import annotations

import json
from pathlib import Path
from typing import Any


class OverlayError(ValueError):
    pass


def load_overlay(path: str | None) -> dict[str, Any]:
    if not path:
        return {}
    with Path(path).open("r", encoding="utf-8") as fh:
        overlay = json.load(fh)
    validate_overlay(overlay)
    return overlay


def validate_overlay(overlay: dict[str, Any]) -> None:
    if not isinstance(overlay, dict):
        raise OverlayError("Overlay must be a JSON object")
    for section in ("record", "sources", "mirrors", "software"):
        if section in overlay and not isinstance(overlay[section], dict):
            raise OverlayError(f"Overlay section {section!r} must be an object")


def source_for_role(overlay: dict[str, Any], file_role: str) -> dict[str, Any] | None:
    sources = overlay.get("sources", {})
    if not isinstance(sources, dict):
        return None
    value = sources.get(file_role)
    return value if isinstance(value, dict) else None


def mirror_for_role(overlay: dict[str, Any], file_role: str) -> dict[str, Any] | None:
    mirrors = overlay.get("mirrors", {})
    if not isinstance(mirrors, dict):
        return None
    value = mirrors.get(file_role)
    return value if isinstance(value, dict) else None


def software_overlay(overlay: dict[str, Any], key: str = "dig_gene_set_extractors") -> dict[str, Any] | None:
    software = overlay.get("software", {})
    if not isinstance(software, dict):
        return None
    value = software.get(key)
    return value if isinstance(value, dict) else None

