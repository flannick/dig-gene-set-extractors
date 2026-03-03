from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

_BUNDLE_RESOURCE_PATHS: dict[str, tuple[str, str | None]] = {
    "ccre_ubiquity_hg19": ("atac/ccre/ccre_ubiquity_hg19.tsv.gz", "hg19"),
    "ccre_ubiquity_hg38": ("atac/ccre/ccre_ubiquity_hg38.tsv.gz", "hg38"),
    "ccre_ubiquity_mm10": ("atac/ccre/ccre_ubiquity_mm10.tsv.gz", "mm10"),
    "encode_ccre_hg19": ("atac/ccre/encode_ccre_hg19.bed.gz", "hg19"),
    "encode_ccre_hg38": ("atac/ccre/encode_ccre_hg38.bed.gz", "hg38"),
    "encode_ccre_mm10": ("atac/ccre/encode_ccre_mm10.bed.gz", "mm10"),
    "atac_reference_profiles_hg19": ("atac/atlas/atac_reference_profiles_hg19.tsv.gz", "hg19"),
    "atac_reference_profiles_hg38": ("atac/atlas/atac_reference_profiles_hg38.tsv.gz", "hg38"),
    "atac_reference_profiles_mm10": ("atac/atlas/atac_reference_profiles_mm10.tsv.gz", "mm10"),
}

_GENOME_BUILD_ALIASES: dict[str, str] = {
    "hg19": "hg19",
    "grch37": "hg19",
    "hg38": "hg38",
    "grch38": "hg38",
    "mm10": "mm10",
    "grcm38": "mm10",
}


def _normalize_genome_build(genome_build: str | None) -> str:
    value = str(genome_build or "").strip().lower()
    return _GENOME_BUILD_ALIASES.get(value, value)


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_json(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object at {path}")
    return payload


def _bundle_checksums(bundle_root: Path) -> dict[str, str]:
    manifest_path = bundle_root / "manifest.json"
    if not manifest_path.exists():
        return {}
    payload = _load_json(manifest_path)
    files = payload.get("files", [])
    if not isinstance(files, list):
        return {}
    out: dict[str, str] = {}
    for row in files:
        if not isinstance(row, dict):
            continue
        rel = str(row.get("path", "")).strip()
        sha = str(row.get("sha256", "")).strip().lower()
        if rel and sha:
            out[rel] = sha
    return out


def _resolve_bundle_file(bundle_root: Path, resource_id: str) -> tuple[Path | None, str | None]:
    rel_and_build = _BUNDLE_RESOURCE_PATHS.get(resource_id)
    if rel_and_build is None:
        return None, None
    rel, genome_build = rel_and_build
    path = bundle_root / rel
    if not path.exists():
        return None, genome_build
    return path, genome_build


def _resource_genome_build(resource_row: dict[str, object], resource_id: str) -> str:
    mapped = _BUNDLE_RESOURCE_PATHS.get(resource_id)
    if mapped is not None:
        _rel, mapped_build = mapped
        return _normalize_genome_build(mapped_build)
    return _normalize_genome_build(str(resource_row.get("genome_build", "")))


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    default_base_manifest = repo_root / "src" / "geneset_extractors" / "resources" / "manifest.json"
    ap = argparse.ArgumentParser(
        description=(
            "Create a local resource manifest from an extracted ATAC reference bundle. "
            "Default layout is direct mode (no fetch step required)."
        )
    )
    ap.add_argument("--bundle-root", required=True, help="Path to extracted bundle directory (contains atac/ and manifest.json)")
    ap.add_argument(
        "--base-manifest",
        default=str(default_base_manifest),
        help="Base resource manifest template to populate",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output manifest path",
    )
    ap.add_argument(
        "--layout",
        choices=["direct", "cache"],
        default="direct",
        help=(
            "direct: set filename to bundle-relative paths and read files directly using --resources_dir <bundle-root> "
            "(no resources fetch step). "
            "cache: keep cache-style filename entries and set url to absolute bundle paths for resources fetch."
        ),
    )
    ap.add_argument(
        "--genome-build",
        "--genome_build",
        dest="genome_build",
        choices=["all", "hg19", "hg38", "mm10"],
        default="all",
        help="Filter resources/presets to a specific genome build (recommended for split bundles)",
    )
    args = ap.parse_args()

    bundle_root = Path(args.bundle_root).expanduser().resolve()
    base_manifest = Path(args.base_manifest).expanduser().resolve()
    out_path = Path(args.out).expanduser().resolve()

    if not bundle_root.exists():
        raise FileNotFoundError(f"Bundle root not found: {bundle_root}")
    if not base_manifest.exists():
        raise FileNotFoundError(f"Base manifest not found: {base_manifest}")

    target_build = _normalize_genome_build(args.genome_build)
    if target_build == "all":
        target_build = ""

    payload = _load_json(base_manifest)
    resources_in = payload.get("resources", [])
    if not isinstance(resources_in, list):
        raise ValueError("Base manifest must contain a resources list")
    checksums = _bundle_checksums(bundle_root)

    resources_out: list[dict[str, object]] = []
    kept_ids: set[str] = set()
    missing: list[str] = []
    resolved: list[str] = []

    for row in resources_in:
        if not isinstance(row, dict):
            continue
        rid = str(row.get("id", "")).strip()
        if not rid:
            continue
        row_build = _resource_genome_build(row, rid)
        if target_build and row_build and row_build != target_build:
            continue
        path, gb = _resolve_bundle_file(bundle_root, rid)
        rel: str | None = None
        if path is not None:
            rel = str(path.relative_to(bundle_root))
        else:
            mapping_path, _mapping_gb = _BUNDLE_RESOURCE_PATHS.get(rid, (None, None))
            rel = mapping_path

        if args.layout == "direct" and rel:
            row["filename"] = rel

        if path is None:
            if gb is not None:
                missing.append(rid)
            row["url"] = "" if args.layout == "direct" else str(row.get("url", ""))
            row["sha256"] = ""
        else:
            if args.layout == "direct":
                row["url"] = ""
            else:
                row["url"] = str(path)
            row["sha256"] = checksums.get(rel, _sha256_file(path))
            resolved.append(rid)

        resources_out.append(row)
        kept_ids.add(rid)

    payload["resources"] = resources_out

    presets = payload.get("presets", {})
    if not isinstance(presets, dict):
        raise ValueError("Base manifest must contain a presets object")
    filtered_presets: dict[str, list[str]] = {}
    for key, value in presets.items():
        if not isinstance(value, list):
            continue
        ids = [str(v).strip() for v in value if str(v).strip() and str(v).strip() in kept_ids]
        if ids:
            filtered_presets[str(key)] = ids
    if target_build:
        per_build_default = f"atac_default_optional_{target_build}"
        if per_build_default in filtered_presets:
            filtered_presets["atac_default_optional"] = list(filtered_presets[per_build_default])
    payload["presets"] = filtered_presets

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"Wrote local resource manifest: {out_path}")
    print(f"Layout: {args.layout}")
    print(f"Genome build filter: {target_build or 'all'}")
    if args.layout == "direct":
        print(f"Use this at runtime with --resources_manifest {out_path} --resources_dir {bundle_root}")
    else:
        print("Use this with `geneset-extractors resources fetch` to stage files into cache.")
    if resolved:
        print("Resolved resource ids:")
        for rid in sorted(resolved):
            print(f"  - {rid}")
    if missing:
        print("Missing bundle files for resource ids (left manual):")
        for rid in sorted(missing):
            print(f"  - {rid}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
