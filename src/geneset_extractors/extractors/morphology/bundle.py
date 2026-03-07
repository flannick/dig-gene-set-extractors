from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
import sys

from geneset_extractors.resource_manager import default_resources_dir, load_manifest, resource_metadata_record


@dataclass
class MorphologyBundleContext:
    manifest_label: str | None = None
    resources: dict[str, dict[str, object]] = field(default_factory=dict)
    resources_dir: Path | None = None
    used: list[dict[str, object]] = field(default_factory=list)
    missing: list[dict[str, object]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def load_bundle_context(resources_manifest: str | None, resources_dir: str | None) -> MorphologyBundleContext:
    manifest_label, resources, _presets, warnings = load_manifest(resources_manifest)
    root = Path(resources_dir).expanduser() if resources_dir else default_resources_dir()
    for warning in warnings:
        print(f"warning: {warning}", file=sys.stderr)
    return MorphologyBundleContext(
        manifest_label=manifest_label,
        resources=resources,
        resources_dir=root,
        warnings=list(warnings),
    )


def resolve_bundle_manifest(
    *,
    ctx: MorphologyBundleContext,
    bundle_id: str,
    resource_policy: str,
) -> tuple[Path | None, dict[str, object] | None]:
    rid = str(bundle_id).strip()
    if ctx.resources_dir is not None:
        local_path = ctx.resources_dir / f"{rid}.bundle.json"
        if local_path.exists():
            ctx.used.append(
                {
                    "id": rid,
                    "method": "reference_bundle_local",
                    "path": str(local_path),
                    "provider": "local_resources_dir",
                    "stable_id": rid,
                    "version": "",
                    "sha256": "",
                    "license": "",
                }
            )
            payload = json.loads(local_path.read_text(encoding="utf-8"))
            payload["_bundle_resolution"] = "local_resources_dir"
            return local_path, payload
    if rid not in ctx.resources:
        raise ValueError(f"Unknown morphology bundle resource id: {rid}")
    entry = ctx.resources[rid]
    path = ctx.resources_dir / str(entry.get("filename", ""))
    if not path.exists():
        ctx.missing.append({"resource_id": rid, "expected_path": str(path), "role_label": "reference_bundle"})
        message = f"Missing morphology reference bundle manifest: {path}. Provide --resources_dir with bundle files or explicit reference files."
        if resource_policy == "fail":
            raise FileNotFoundError(message)
        print(f"warning: {message}", file=sys.stderr)
        return None, None
    ctx.used.append(resource_metadata_record(resource_id=rid, entry=entry, local_path=path, method="reference_bundle"))
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["_bundle_resolution"] = "resource_manifest"
    return path, payload


def bundle_resources_info(ctx: MorphologyBundleContext | None) -> dict[str, object] | None:
    if ctx is None:
        return None
    if not ctx.used and not ctx.missing and not ctx.warnings:
        return None
    return {
        "manifest": ctx.manifest_label,
        "resources_dir": str(ctx.resources_dir) if ctx.resources_dir else None,
        "used": list(ctx.used),
        "missing": list(ctx.missing),
        "warnings": list(ctx.warnings),
    }
