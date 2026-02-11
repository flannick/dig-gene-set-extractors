from __future__ import annotations

import json
from pathlib import Path

from omics2geneset.converters import atac_bulk, atac_sc_10x


CONVERTERS = {
    "atac_bulk": atac_bulk,
    "atac_sc_10x": atac_sc_10x,
}


def _schemas_dir() -> Path:
    return Path(__file__).parent / "schemas"


def _specs_dir() -> Path:
    return Path(__file__).parent / "converters" / "specs"


def _validate_converter_spec(payload: dict[str, object]) -> None:
    required = ["name", "description", "inputs", "parameters", "outputs"]
    missing = [k for k in required if k not in payload]
    if missing:
        raise ValueError(f"converter spec missing required fields: {', '.join(missing)}")


def list_converters() -> list[tuple[str, str]]:
    out = []
    for name in sorted(CONVERTERS):
        spec = get_converter_spec(name)
        out.append((name, str(spec.get("description", ""))))
    return out


def get_converter(name: str):
    if name not in CONVERTERS:
        raise KeyError(f"Unknown converter: {name}")
    return CONVERTERS[name]


def get_converter_spec(name: str) -> dict[str, object]:
    p = _specs_dir() / f"{name}.json"
    if not p.exists():
        raise FileNotFoundError(f"Missing converter spec: {p}")
    payload = json.loads(p.read_text(encoding="utf-8"))
    _validate_converter_spec(payload)
    return payload
