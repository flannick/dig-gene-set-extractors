from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import sys

from geneset_extractors.resource_manager import default_resources_dir, load_manifest, resource_metadata_record


@dataclass
class ResourceContext:
    manifest_label: str | None = None
    resources: dict[str, dict[str, object]] = field(default_factory=dict)
    resources_dir: Path | None = None
    warnings: list[str] = field(default_factory=list)
    used: list[dict[str, object]] = field(default_factory=list)
    missing: list[dict[str, object]] = field(default_factory=list)



def load_resource_context(resources_manifest: str | None, resources_dir: str | None) -> ResourceContext:
    root = Path(resources_dir).expanduser() if resources_dir else default_resources_dir()
    manifest_path = resources_manifest
    if manifest_path is None:
        local_manifest = root / "local_resources_manifest.json"
        if local_manifest.exists():
            manifest_path = str(local_manifest)

    manifest_label, resources, _presets, warnings = load_manifest(manifest_path)
    for warning in warnings:
        print(f"warning: {warning}", file=sys.stderr)
    return ResourceContext(
        manifest_label=manifest_label,
        resources=resources,
        resources_dir=root,
        warnings=list(warnings),
        used=[],
        missing=[],
    )



def resolve_resource_path(*, ctx: ResourceContext, resource_id: str | None, resource_policy: str, role_label: str, enablement_hint: str) -> Path | None:
    rid = "" if resource_id is None else str(resource_id).strip()
    if not rid:
        return None
    if rid not in ctx.resources:
        raise ValueError(f"Unknown {role_label} resource id: {rid}")
    entry = ctx.resources[rid]
    filename = str(entry.get("filename", ""))
    candidates = [ctx.resources_dir / filename, ctx.resources_dir / "bundle" / filename]
    path = next((candidate for candidate in candidates if candidate.exists()), None)
    if path is not None:
        ctx.used.append(resource_metadata_record(resource_id=rid, entry=entry, local_path=path, method=role_label))
        return path

    message = (
        f"Missing {role_label} resource file. Tried: "
        + ", ".join(str(candidate) for candidate in candidates)
        + f". {enablement_hint}"
    )
    ctx.missing.append({"resource_id": rid, "expected_path": str(candidates[0]), "role_label": role_label})
    if resource_policy == "fail":
        raise FileNotFoundError(message)
    print(f"warning: {message}", file=sys.stderr)
    return None



def build_resources_info(ctx: ResourceContext) -> dict[str, object] | None:
    if not ctx.used and not ctx.missing and not ctx.warnings:
        return None
    return {
        "manifest": ctx.manifest_label,
        "resources_dir": str(ctx.resources_dir) if ctx.resources_dir else None,
        "used": list(ctx.used),
        "missing": list(ctx.missing),
        "warnings": list(ctx.warnings),
    }
