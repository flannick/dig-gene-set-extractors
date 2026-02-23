from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import tarfile
import tempfile
from typing import Iterable


_COMMON_FILES = [
    "LICENSES/THIRD_PARTY_NOTICES.txt",
    "atac/ccre/core_collection_list.txt.gz",
    "atac/ccre/core_collection_class_matrix.txt.gz",
]

_BUILD_FILES = {
    "hg19": [
        "atac/ccre/encode_ccre_hg19.bed.gz",
        "atac/ccre/ccre_ubiquity_hg19.tsv.gz",
        "atac/gtf/gencode.v19.annotation.gtf.gz",
        "atac/atlas/atac_reference_profiles_hg19.tsv.gz",
    ],
    "hg38": [
        "atac/ccre/encode_ccre_hg38.bed.gz",
        "atac/ccre/ccre_ubiquity_hg38.tsv.gz",
        "atac/gtf/gencode.v49.basic.annotation.gtf.gz",
        "atac/atlas/atac_reference_profiles_hg38.tsv.gz",
    ],
}


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_manifest(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object at {path}")
    return payload


def _copy_files(bundle_src: Path, bundle_dst: Path, rel_paths: Iterable[str]) -> None:
    for rel in rel_paths:
        src = bundle_src / rel
        if not src.exists():
            raise FileNotFoundError(f"Missing source bundle file: {src}")
        dst = bundle_dst / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def _source_url_index(source_manifest: dict[str, object]) -> dict[str, object]:
    out: dict[str, object] = {}
    files = source_manifest.get("files", [])
    if not isinstance(files, list):
        return out
    for row in files:
        if not isinstance(row, dict):
            continue
        rel = str(row.get("path", "")).strip()
        if not rel:
            continue
        out[rel] = row.get("source_url")
    return out


def _build_manifest(
    bundle_name: str,
    version: str,
    bundle_dir: Path,
    rel_paths: list[str],
    source_manifest: dict[str, object],
) -> None:
    source_url_by_path = _source_url_index(source_manifest)
    file_rows: list[dict[str, object]] = []
    source_url_strings: set[str] = set()

    for rel in sorted(rel_paths + ["README.txt"]):
        p = bundle_dir / rel
        row = {
            "path": rel,
            "size_bytes": p.stat().st_size,
            "sha256": _sha256_file(p),
            "source_url": source_url_by_path.get(rel),
        }
        file_rows.append(row)
        src = row["source_url"]
        if isinstance(src, str) and src:
            source_url_strings.add(src)
        if isinstance(src, list):
            for x in src:
                if isinstance(x, str) and x:
                    source_url_strings.add(x)

    raw_inputs_out: list[dict[str, object]] = []
    raw_inputs = source_manifest.get("raw_inputs", [])
    if isinstance(raw_inputs, list):
        for row in raw_inputs:
            if not isinstance(row, dict):
                continue
            src = row.get("source_url")
            src_values = [src] if isinstance(src, str) else src if isinstance(src, list) else []
            if any(isinstance(x, str) and x in source_url_strings for x in src_values):
                raw_inputs_out.append(dict(row))

    payload = {
        "bundle_version": version,
        "bundle_name": bundle_name,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "bundle_root": "bundle",
        "raw_inputs": raw_inputs_out,
        "files": file_rows,
        "qc": {
            "source_manifest_bundle_name": source_manifest.get("bundle_name"),
            "source_manifest_bundle_version": source_manifest.get("bundle_version"),
        },
    }
    (bundle_dir / "manifest.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _write_readme(bundle_dir: Path, build: str, version: str) -> None:
    text = (
        f"DIG Gene Set Extractors ATAC Reference Bundle ({build}) v{version}\n\n"
        "This is a build-specific bundle for human ATAC reference resources.\n\n"
        "Included resources:\n"
        f"- encode_ccre_{build}.bed.gz\n"
        f"- ccre_ubiquity_{build}.tsv.gz\n"
        f"- atac_reference_profiles_{build}.tsv.gz\n"
        "- core_collection_list.txt.gz\n"
        "- core_collection_class_matrix.txt.gz\n\n"
        "Use this bundle with scripts/make_local_resources_manifest.py and set\n"
        f"--genome-build {build} when generating the local resource manifest.\n"
    )
    (bundle_dir / "README.txt").write_text(text, encoding="utf-8")


def _create_tarball(bundle_root: Path, out_tarball: Path) -> None:
    if out_tarball.exists():
        out_tarball.unlink()
    bundle_dir = bundle_root / "bundle"
    with tarfile.open(out_tarball, "w:gz") as tf:
        tf.add(bundle_dir, arcname="bundle")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Split a combined human ATAC reference bundle into hg19/hg38 build-specific tarballs."
    )
    ap.add_argument(
        "--source-bundle",
        default="refdata_build/bundle",
        help="Path to source bundle directory (contains manifest.json and atac/)",
    )
    ap.add_argument(
        "--out-dir",
        default="refdata_build",
        help="Output directory for build-specific tarballs",
    )
    ap.add_argument(
        "--builds",
        default="hg19,hg38",
        help="Comma-separated builds to emit (subset of hg19,hg38)",
    )
    ap.add_argument(
        "--version",
        default="",
        help="Bundle version override (defaults to source manifest bundle_version)",
    )
    args = ap.parse_args()

    source_bundle = Path(args.source_bundle).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    source_manifest_path = source_bundle / "manifest.json"
    if not source_manifest_path.exists():
        raise FileNotFoundError(f"Missing source manifest: {source_manifest_path}")
    source_manifest = _load_manifest(source_manifest_path)
    version = str(args.version).strip() or str(source_manifest.get("bundle_version", "")).strip()
    if not version:
        raise ValueError("Could not resolve bundle version (pass --version)")

    requested_builds = [b.strip().lower() for b in str(args.builds).split(",") if b.strip()]
    for build in requested_builds:
        if build not in _BUILD_FILES:
            raise ValueError(f"Unsupported build for split: {build}")

    checksums: list[tuple[str, str]] = []
    for build in requested_builds:
        rel_paths = list(_COMMON_FILES) + list(_BUILD_FILES[build])
        bundle_name = f"dig-gene-set-extractors-atac-refdata-{build}-v{version}"
        out_tarball = out_dir / f"{bundle_name}.tar.gz"

        with tempfile.TemporaryDirectory(prefix=f"split-{build}-", dir=str(out_dir)) as tmp:
            tmp_root = Path(tmp)
            bundle_dir = tmp_root / "bundle"
            bundle_dir.mkdir(parents=True, exist_ok=True)
            _copy_files(source_bundle, bundle_dir, rel_paths)
            _write_readme(bundle_dir, build, version)
            _build_manifest(bundle_name, version, bundle_dir, rel_paths, source_manifest)
            _create_tarball(tmp_root, out_tarball)

        checksums.append((_sha256_file(out_tarball), out_tarball.name))
        print(f"Wrote {out_tarball}")

    sha_path = out_dir / "SHA256SUMS.human_split.txt"
    with sha_path.open("w", encoding="utf-8") as fh:
        for sha, name in checksums:
            fh.write(f"{sha}  {name}\n")
    print(f"Wrote {sha_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
