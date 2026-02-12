from __future__ import annotations

import csv
import json
from pathlib import Path


def _validate_required_fallback(payload: object, schema: dict[str, object], path: str = "$") -> None:
    if not isinstance(payload, dict):
        raise ValueError(f"{path}: expected object")
    required = schema.get("required", [])
    if isinstance(required, list):
        for key in required:
            if key not in payload:
                raise ValueError(f"{path}: missing required key '{key}'")
    props = schema.get("properties", {})
    if not isinstance(props, dict):
        return
    for key, sub_schema in props.items():
        if key in payload and isinstance(sub_schema, dict) and sub_schema.get("type") == "object":
            _validate_required_fallback(payload[key], sub_schema, f"{path}.{key}")


def validate_metadata_schema(meta_path: Path, schema_path: Path) -> None:
    payload = json.loads(meta_path.read_text(encoding="utf-8"))
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    try:
        import jsonschema  # type: ignore
    except ModuleNotFoundError:
        # Compatibility fallback when jsonschema is not yet installed.
        _validate_required_fallback(payload, schema)
    else:
        jsonschema.validate(payload, schema)


def validate_geneset_tsv(geneset_path: Path) -> None:
    with geneset_path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"gene_id", "weight"}
        if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
            raise ValueError("geneset.tsv missing required columns gene_id and weight")
        seen: set[str] = set()
        for row in reader:
            gid = row["gene_id"]
            if gid in seen:
                raise ValueError(f"duplicate gene_id in geneset.tsv: {gid}")
            seen.add(gid)
            float(row["weight"])


def _resolve_manifest_path(root: Path, manifest_value: str) -> Path:
    candidate = Path(manifest_value)
    if candidate.is_absolute():
        return candidate
    joined = root / candidate
    if joined.exists():
        return joined
    return candidate


def _validate_grouped_output_dir(out_dir: Path, schema_path: Path) -> dict[str, object]:
    manifest = out_dir / "manifest.tsv"
    if not manifest.exists():
        raise FileNotFoundError("output dir must contain geneset.tsv and geneset.meta.json, or manifest.tsv for grouped outputs")
    with manifest.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or "path" not in reader.fieldnames:
            raise ValueError("manifest.tsv must contain a path column")
        rows = list(reader)
    if not rows:
        raise ValueError("manifest.tsv contains no group rows")
    failures: list[str] = []
    for row in rows:
        group_path = _resolve_manifest_path(out_dir, str(row["path"]))
        try:
            _validate_single_output_dir(group_path, schema_path)
        except Exception as exc:
            failures.append(f"{group_path}: {exc}")
    if failures:
        raise ValueError("grouped validation failed: " + "; ".join(failures))
    return {"mode": "grouped", "n_groups": len(rows)}


def _validate_single_output_dir(out: Path, schema_path: Path) -> None:
    geneset = out / "geneset.tsv"
    meta = out / "geneset.meta.json"
    if not geneset.exists() or not meta.exists():
        raise FileNotFoundError("output dir must contain geneset.tsv and geneset.meta.json")
    validate_geneset_tsv(geneset)
    validate_metadata_schema(meta, schema_path)


def validate_output_dir(out_dir: str | Path, schema_path: str | Path) -> dict[str, object]:
    out = Path(out_dir)
    schema = Path(schema_path)
    geneset = out / "geneset.tsv"
    meta = out / "geneset.meta.json"
    if geneset.exists() and meta.exists():
        _validate_single_output_dir(out, schema)
        return {"mode": "single", "n_groups": 1}
    return _validate_grouped_output_dir(out, schema)
