from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
from urllib.parse import urlparse
from urllib.request import url2pathname, urlopen


_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")

_MANIFEST_SCHEMA: dict[str, object] = {
    "type": "object",
    "required": ["resources", "presets"],
    "properties": {
        "resources": {
            "type": "array",
            "items": {
                "type": "object",
                "required": ["id", "filename"],
                "properties": {
                    "id": {"type": "string"},
                    "description": {"type": "string"},
                    "provider": {"type": "string"},
                    "stable_id": {"type": "string"},
                    "version": {"type": "string"},
                    "genome_build": {"type": "string"},
                    "filename": {"type": "string"},
                    "url": {"type": "string"},
                    "sha256": {"type": "string"},
                    "license": {"type": "string"},
                },
                "additionalProperties": True,
            },
        },
        "presets": {
            "type": "object",
            "additionalProperties": {
                "type": "array",
                "items": {"type": "string"},
            },
        },
    },
    "additionalProperties": True,
}


def default_resources_dir() -> Path:
    env = os.getenv("OMICS2GENESET_RESOURCES_DIR")
    if env:
        return Path(env).expanduser()
    return Path.home() / ".cache" / "omics2geneset" / "resources"


def bundled_manifest_path() -> Path:
    return Path(__file__).parent / "resources" / "manifest.json"


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _read_manifest_payload(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:  # pragma: no cover
        raise ValueError(f"Could not parse manifest JSON at {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"Manifest at {path} must be a JSON object")
    return payload


def _validate_manifest_schema(payload: dict[str, object], manifest_label: str) -> None:
    try:
        import jsonschema  # type: ignore
    except ModuleNotFoundError:
        if "resources" not in payload or "presets" not in payload:
            raise ValueError(f"Manifest {manifest_label} must include resources and presets")
        return
    jsonschema.validate(payload, _MANIFEST_SCHEMA)


def _merge_manifest_payloads(base: dict[str, object], overlay: dict[str, object]) -> dict[str, object]:
    base_resources = base.get("resources", [])
    overlay_resources = overlay.get("resources", [])
    if not isinstance(base_resources, list) or not isinstance(overlay_resources, list):
        raise ValueError("resource manifest: resources must be a list")

    merged_resources: list[dict[str, object]] = []
    index_by_id: dict[str, int] = {}
    for entry in base_resources:
        if not isinstance(entry, dict):
            raise ValueError("resource manifest: each resource entry must be an object")
        rid = str(entry.get("id", "")).strip()
        if not rid:
            raise ValueError("resource manifest: each resource entry must include non-empty id")
        index_by_id[rid] = len(merged_resources)
        merged_resources.append(dict(entry))

    for entry in overlay_resources:
        if not isinstance(entry, dict):
            raise ValueError("resource manifest: each resource entry must be an object")
        rid = str(entry.get("id", "")).strip()
        if not rid:
            raise ValueError("resource manifest: each resource entry must include non-empty id")
        normalized = dict(entry)
        if rid in index_by_id:
            merged_resources[index_by_id[rid]] = normalized
        else:
            index_by_id[rid] = len(merged_resources)
            merged_resources.append(normalized)

    base_presets = base.get("presets", {})
    overlay_presets = overlay.get("presets", {})
    if not isinstance(base_presets, dict) or not isinstance(overlay_presets, dict):
        raise ValueError("resource manifest: presets must be an object")
    merged_presets: dict[str, list[str]] = {}
    for key, value in base_presets.items():
        if isinstance(value, list):
            merged_presets[str(key)] = [str(v) for v in value]
    for key, value in overlay_presets.items():
        if not isinstance(value, list):
            raise ValueError(f"resource manifest: preset {key} must be a list")
        merged_presets[str(key)] = [str(v) for v in value]

    return {
        "resources": merged_resources,
        "presets": merged_presets,
    }


def _extract_manifest_data(
    payload: dict[str, object],
) -> tuple[dict[str, dict[str, object]], dict[str, list[str]], list[str]]:
    resources_raw = payload.get("resources", [])
    if not isinstance(resources_raw, list):
        raise ValueError("resource manifest: resources must be a list")

    resources: dict[str, dict[str, object]] = {}
    warnings: list[str] = []
    for entry in resources_raw:
        if not isinstance(entry, dict):
            raise ValueError("resource manifest: each resource entry must be an object")
        rid = str(entry.get("id", "")).strip()
        if not rid:
            raise ValueError("resource manifest: each resource entry must include non-empty id")
        if rid in resources:
            raise ValueError(f"resource manifest: duplicate resource id: {rid}")
        filename = str(entry.get("filename", "")).strip()
        if not filename:
            raise ValueError(f"resource manifest: resource {rid} missing filename")

        normalized = dict(entry)
        normalized["id"] = rid
        normalized["filename"] = filename
        normalized["url"] = str(normalized.get("url", "")).strip()
        normalized["sha256"] = str(normalized.get("sha256", "")).strip().lower()
        normalized["stable_id"] = str(normalized.get("stable_id", "")).strip()
        normalized["version"] = str(normalized.get("version", "")).strip()
        normalized["license"] = str(normalized.get("license", "")).strip()
        normalized["provider"] = str(normalized.get("provider", "")).strip()
        normalized["description"] = str(normalized.get("description", "")).strip()
        normalized["genome_build"] = str(normalized.get("genome_build", "")).strip()
        resources[rid] = normalized

        url = str(normalized.get("url", "")).strip()
        sha = str(normalized.get("sha256", "")).strip()
        if sha and not _SHA256_RE.match(sha):
            raise ValueError(f"resource manifest: resource {rid} has invalid sha256 (must be 64 lowercase hex chars)")
        if url and not sha:
            warnings.append(f"resource {rid} has url but no sha256 (download is unverified)")
        if url and not str(normalized.get("stable_id", "")).strip():
            warnings.append(f"resource {rid} has url but no stable_id")
        if url and not str(normalized.get("version", "")).strip():
            warnings.append(f"resource {rid} has url but no version")

    presets_raw = payload.get("presets", {})
    if not isinstance(presets_raw, dict):
        raise ValueError("resource manifest: presets must be an object")
    presets: dict[str, list[str]] = {}
    for k, v in presets_raw.items():
        key = str(k).strip()
        if not key:
            continue
        if not isinstance(v, list):
            raise ValueError(f"resource manifest: preset {key} must be a list")
        out_ids: list[str] = []
        seen: set[str] = set()
        for item in v:
            rid = str(item).strip()
            if not rid:
                continue
            if rid not in resources:
                raise ValueError(f"resource manifest: preset {key} references unknown resource id: {rid}")
            if rid in seen:
                continue
            out_ids.append(rid)
            seen.add(rid)
        presets[key] = out_ids
    return resources, presets, warnings


def load_manifest(
    manifest_path: str | Path | None = None,
    merge_with_bundled: bool = True,
) -> tuple[str, dict[str, dict[str, object]], dict[str, list[str]], list[str]]:
    bundled = bundled_manifest_path()
    if manifest_path is None:
        payload = _read_manifest_payload(bundled)
        manifest_label = str(bundled)
    else:
        user_path = Path(manifest_path).expanduser()
        if merge_with_bundled and user_path.resolve() != bundled.resolve():
            payload = _merge_manifest_payloads(
                _read_manifest_payload(bundled),
                _read_manifest_payload(user_path),
            )
            manifest_label = f"{bundled}+{user_path}"
        else:
            payload = _read_manifest_payload(user_path)
            manifest_label = str(user_path)

    _validate_manifest_schema(payload, manifest_label)
    resources, presets, warnings = _extract_manifest_data(payload)
    return manifest_label, resources, presets, warnings


def resolve_requested_resource_ids(
    resources: dict[str, dict[str, object]],
    presets: dict[str, list[str]],
    resource_ids: list[str],
    preset: str | None,
) -> list[str]:
    requested: list[str] = []
    if preset is not None:
        p = str(preset).strip()
        if not p:
            raise ValueError("preset must be non-empty when provided")
        if p not in presets:
            raise ValueError(f"Unknown resource preset: {p}")
        requested.extend(presets[p])
    requested.extend([str(r).strip() for r in resource_ids if str(r).strip()])

    out: list[str] = []
    seen: set[str] = set()
    for rid in requested:
        if rid not in resources:
            raise ValueError(f"Unknown resource id: {rid}")
        if rid in seen:
            continue
        out.append(rid)
        seen.add(rid)
    if not out:
        raise ValueError("No resource ids selected. Provide ids and/or --preset.")
    return out


def resource_status_rows(
    resources: dict[str, dict[str, object]],
    resources_dir: str | Path | None = None,
    verify: bool = False,
) -> list[dict[str, object]]:
    root = default_resources_dir() if resources_dir is None else Path(resources_dir).expanduser()
    rows: list[dict[str, object]] = []
    for rid in sorted(resources):
        entry = resources[rid]
        filename = str(entry.get("filename", ""))
        path = root / filename
        expected_sha = str(entry.get("sha256", "")).strip().lower()
        url = str(entry.get("url", "")).strip()
        availability = "download" if url else "manual"
        exists = path.exists()
        size_bytes = path.stat().st_size if exists else None

        should_hash = bool(exists and (verify or expected_sha))
        actual_sha = sha256_path(path) if should_hash else None

        if not exists:
            status = "manual_missing" if availability == "manual" else "missing"
        elif expected_sha:
            status = "ok" if actual_sha == expected_sha else "checksum_mismatch"
        elif verify:
            status = "ok_unverified"
        else:
            status = "ok_fast"

        rows.append(
            {
                "id": rid,
                "status": status,
                "path": str(path),
                "exists": exists,
                "size_bytes": size_bytes,
                "availability": availability,
                "expected_sha256": expected_sha or None,
                "actual_sha256": actual_sha,
                "url": url or None,
            }
        )
    return rows


def _download_to_path(url: str, out_path: Path) -> None:
    parsed = urlparse(url)
    if parsed.scheme in {"", "file"}:
        if parsed.scheme == "file":
            src = Path(url2pathname(parsed.path))
        else:
            src = Path(url)
        if not src.exists():
            raise FileNotFoundError(f"Resource source does not exist: {src}")
        shutil.copyfile(src, out_path)
        return

    with urlopen(url) as response, out_path.open("wb") as out_fh:
        shutil.copyfileobj(response, out_fh)


def fetch_resources(
    resources: dict[str, dict[str, object]],
    resource_ids: list[str],
    resources_dir: str | Path | None = None,
    overwrite: bool = False,
    skip_missing_urls: bool = True,
) -> list[dict[str, object]]:
    root = default_resources_dir() if resources_dir is None else Path(resources_dir).expanduser()
    root.mkdir(parents=True, exist_ok=True)
    out: list[dict[str, object]] = []

    for rid in resource_ids:
        entry = resources[rid]
        filename = str(entry.get("filename", "")).strip()
        url = str(entry.get("url", "")).strip()
        expected_sha = str(entry.get("sha256", "")).strip().lower()
        target = root / filename

        if not url:
            if skip_missing_urls:
                out.append(
                    {
                        "id": rid,
                        "status": "manual",
                        "path": str(target),
                        "reason": "missing_url",
                    }
                )
                continue
            raise ValueError(f"Resource {rid} has no url in manifest; cannot fetch.")

        if target.exists() and not overwrite:
            actual = sha256_path(target) if expected_sha else None
            if expected_sha and actual != expected_sha:
                raise ValueError(
                    f"Existing file checksum mismatch for {rid}: {target}. "
                    "Use --overwrite to replace."
                )
            out.append(
                {
                    "id": rid,
                    "status": "skipped",
                    "path": str(target),
                }
            )
            continue

        with tempfile.NamedTemporaryFile(prefix="omics2geneset_resource_", dir=root, delete=False) as tmp:
            tmp_path = Path(tmp.name)
        try:
            _download_to_path(url, tmp_path)
            actual_sha = sha256_path(tmp_path)
            if expected_sha and actual_sha != expected_sha:
                raise ValueError(
                    f"Checksum mismatch for {rid}: expected {expected_sha} got {actual_sha}"
                )
            target.parent.mkdir(parents=True, exist_ok=True)
            tmp_path.replace(target)
        finally:
            if tmp_path.exists():
                tmp_path.unlink()

        out.append(
            {
                "id": rid,
                "status": "fetched" if expected_sha else "fetched_unverified",
                "path": str(target),
            }
        )
    return out


def describe_resource(
    resources: dict[str, dict[str, object]],
    resource_id: str,
    resources_dir: str | Path | None = None,
    verify: bool = False,
) -> dict[str, object]:
    rid = str(resource_id).strip()
    if rid not in resources:
        raise ValueError(f"Unknown resource id: {rid}")
    status = next(row for row in resource_status_rows(resources, resources_dir=resources_dir, verify=verify) if row["id"] == rid)
    entry = dict(resources[rid])
    entry["id"] = rid
    entry["status"] = status["status"]
    entry["path"] = status["path"]
    entry["exists"] = status["exists"]
    entry["size_bytes"] = status["size_bytes"]
    entry["availability"] = status["availability"]
    entry["actual_sha256"] = status["actual_sha256"]
    return entry


def resource_metadata_record(
    resource_id: str,
    entry: dict[str, object],
    local_path: str | Path,
    method: str,
) -> dict[str, object]:
    return {
        "id": resource_id,
        "method": method,
        "path": str(local_path),
        "provider": str(entry.get("provider", "")),
        "stable_id": str(entry.get("stable_id", "")),
        "version": str(entry.get("version", "")),
        "sha256": str(entry.get("sha256", "")),
        "license": str(entry.get("license", "")),
    }
