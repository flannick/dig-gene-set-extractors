from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import sys

from omics2geneset.resource_manager import default_resources_dir, load_manifest, resource_metadata_record


@dataclass
class ResourceContext:
    manifest_label: str | None = None
    resources: dict[str, dict[str, object]] = field(default_factory=dict)
    resources_dir: Path | None = None
    warnings: list[str] = field(default_factory=list)
    used: list[dict[str, object]] = field(default_factory=list)
    missing: list[dict[str, object]] = field(default_factory=list)


def load_resource_context(
    resources_manifest: str | None,
    resources_dir: str | None,
) -> ResourceContext:
    manifest_label, resources, _presets, warnings = load_manifest(resources_manifest)
    root = Path(resources_dir).expanduser() if resources_dir else default_resources_dir()
    for warning in warnings:
        print(f"warning: {warning}", file=sys.stderr)
    return ResourceContext(
        manifest_label=manifest_label,
        resources=resources,
        resources_dir=root,
        warnings=list(warnings),
        used=[],
    )


def resolve_resource_path(
    *,
    ctx: ResourceContext,
    resource_id: str | None,
    resource_policy: str,
    role_label: str,
    enablement_hint: str,
) -> Path | None:
    rid = "" if resource_id is None else str(resource_id).strip()
    if not rid:
        return None
    if rid not in ctx.resources:
        raise ValueError(f"Unknown {role_label} resource id: {rid}")
    entry = ctx.resources[rid]
    path = ctx.resources_dir / str(entry.get("filename", ""))
    if path.exists():
        ctx.used.append(
            resource_metadata_record(
                resource_id=rid,
                entry=entry,
                local_path=path,
                method=role_label,
            )
        )
        return path

    message = (
        f"Missing {role_label} resource file: {path}. {enablement_hint}"
    )
    ctx.missing.append(
        {
            "resource_id": rid,
            "expected_path": str(path),
            "role_label": role_label,
        }
    )
    if resource_policy == "fail":
        raise FileNotFoundError(message)
    print(f"warning: {message}", file=sys.stderr)
    return None


def build_resources_info(ctx: ResourceContext) -> dict[str, object] | None:
    if not ctx.used and not ctx.warnings and not ctx.missing:
        return None
    return {
        "manifest": ctx.manifest_label,
        "resources_dir": str(ctx.resources_dir) if ctx.resources_dir else None,
        "used": list(ctx.used),
        "missing": list(ctx.missing),
        "warnings": list(ctx.warnings),
    }
