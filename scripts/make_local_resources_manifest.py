from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


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
    mapping: dict[str, tuple[str, str | None]] = {
        "ccre_ubiquity_hg38": ("atac/ccre/ccre_ubiquity_hg38.tsv.gz", "hg38"),
        "ccre_ubiquity_mm10": ("atac/ccre/ccre_ubiquity_mm10.tsv.gz", "mm10"),
        "encode_ccre_hg38": ("atac/ccre/encode_ccre_hg38.bed.gz", "hg38"),
        "encode_ccre_mm10": ("atac/ccre/encode_ccre_mm10.bed.gz", "mm10"),
        "atac_reference_profiles_hg38": ("atac/atlas/atac_reference_profiles_hg38.tsv.gz", "hg38"),
        "atac_reference_profiles_mm10": ("atac/atlas/atac_reference_profiles_mm10.tsv.gz", "mm10"),
    }
    rel_and_build = mapping.get(resource_id)
    if rel_and_build is None:
        return None, None
    rel, genome_build = rel_and_build
    path = bundle_root / rel
    if not path.exists():
        return None, genome_build
    return path, genome_build


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Create a resource manifest that points omics2geneset resources to local bundle files "
            "using Unix filesystem paths (no URI required)."
        )
    )
    ap.add_argument("--bundle-root", required=True, help="Path to extracted bundle directory (contains atac/ and manifest.json)")
    ap.add_argument(
        "--base-manifest",
        default="src/omics2geneset/resources/manifest.json",
        help="Base resource manifest template to populate",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output manifest path",
    )
    args = ap.parse_args()

    bundle_root = Path(args.bundle_root).expanduser().resolve()
    base_manifest = Path(args.base_manifest).expanduser().resolve()
    out_path = Path(args.out).expanduser().resolve()

    if not bundle_root.exists():
        raise FileNotFoundError(f"Bundle root not found: {bundle_root}")
    if not base_manifest.exists():
        raise FileNotFoundError(f"Base manifest not found: {base_manifest}")

    payload = _load_json(base_manifest)
    resources = payload.get("resources", [])
    if not isinstance(resources, list):
        raise ValueError("Base manifest must contain a resources list")
    checksums = _bundle_checksums(bundle_root)

    missing: list[str] = []
    resolved: list[str] = []

    for row in resources:
        if not isinstance(row, dict):
            continue
        rid = str(row.get("id", "")).strip()
        if not rid:
            continue
        path, gb = _resolve_bundle_file(bundle_root, rid)
        if path is None:
            if gb in {"hg38", "mm10"}:
                missing.append(rid)
            row["url"] = ""
            row["sha256"] = ""
            continue
        rel = str(path.relative_to(bundle_root))
        row["url"] = str(path)
        row["sha256"] = checksums.get(rel, _sha256_file(path))
        resolved.append(rid)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"Wrote local resource manifest: {out_path}")
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
