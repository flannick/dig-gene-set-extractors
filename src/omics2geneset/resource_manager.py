from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
from urllib.parse import urlparse
from urllib.request import url2pathname, urlopen


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


def load_manifest(manifest_path: str | Path | None = None) -> tuple[Path, dict[str, dict[str, object]], dict[str, list[str]]]:
    path = Path(manifest_path) if manifest_path is not None else bundled_manifest_path()
    payload = json.loads(path.read_text(encoding="utf-8"))
    resources_raw = payload.get("resources", [])
    if not isinstance(resources_raw, list):
        raise ValueError("resource manifest: resources must be a list")
    resources: dict[str, dict[str, object]] = {}
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
        resources[rid] = entry

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
        ids = [str(x).strip() for x in v if str(x).strip()]
        presets[key] = ids
    return path, resources, presets


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
) -> list[dict[str, object]]:
    root = default_resources_dir() if resources_dir is None else Path(resources_dir).expanduser()
    rows: list[dict[str, object]] = []
    for rid in sorted(resources):
        entry = resources[rid]
        filename = str(entry.get("filename", ""))
        path = root / filename
        expected_sha = str(entry.get("sha256", "")).strip().lower()
        exists = path.exists()
        actual_sha = sha256_path(path) if exists else None
        if not exists:
            status = "missing"
        elif expected_sha and actual_sha != expected_sha:
            status = "checksum_mismatch"
        else:
            status = "ok"
        rows.append(
            {
                "id": rid,
                "status": status,
                "path": str(path),
                "exists": exists,
                "expected_sha256": expected_sha,
                "actual_sha256": actual_sha,
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

        if target.exists() and not overwrite:
            actual = sha256_path(target)
            if expected_sha and actual != expected_sha:
                raise ValueError(
                    f"Existing file checksum mismatch for {rid}: {target}. "
                    "Use --overwrite to replace."
                )
            out.append({"id": rid, "status": "skipped", "path": str(target)})
            continue

        if not url:
            raise ValueError(f"Resource {rid} has no url in manifest; cannot fetch.")

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

        out.append({"id": rid, "status": "fetched", "path": str(target)})
    return out
