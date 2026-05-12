from __future__ import annotations

from pathlib import Path
from typing import Any


def write_reproduce_sh(record: dict[str, Any], out_path: Path) -> None:
    replay = record.get("replay", {})
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        'ROOT="${1:-$PWD/replay}"',
        'mkdir -p "$ROOT"',
        'cd "$ROOT"',
        "",
    ]
    for step in replay.get("steps", []):
        if not isinstance(step, dict):
            continue
        cmd = step.get("replay_command") or step.get("observed_command")
        if not cmd:
            continue
        lines.append(f'echo "Running {step.get("name", step.get("id", "step"))}"')
        lines.append(str(cmd))
        lines.append("")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    out_path.chmod(0o755)


def _iter_replay_files(record: dict[str, Any]) -> list[dict[str, Any]]:
    replay = record.get("replay", {})
    files: list[dict[str, Any]] = []
    for item in replay.get("inputs", []):
        if isinstance(item, dict):
            files.append(item)
    for step in replay.get("steps", []):
        if not isinstance(step, dict):
            continue
        for item in step.get("generated_files", []):
            if isinstance(item, dict):
                files.append(item)
    return files


def write_checksums(record: dict[str, Any], out_path: Path, role_filter: set[str] | None = None) -> None:
    rows: list[str] = []
    for item in _iter_replay_files(record):
        role = str(item.get("role", ""))
        if role_filter is not None and role not in role_filter:
            continue
        rel = item.get("path")
        if not rel:
            continue
        for checksum in item.get("checksums", []):
            if not isinstance(checksum, dict):
                continue
            rows.append(
                "\t".join(
                    [
                        str(rel),
                        role,
                        str(checksum.get("algorithm", "")),
                        str(checksum.get("encoding", "")),
                        str(checksum.get("value", "")),
                    ]
                )
            )
    out_path.write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")


def write_input_downloads(record: dict[str, Any], out_path: Path) -> None:
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", ""]
    replay = record.get("replay", {})
    for item in replay.get("inputs", []):
        if not isinstance(item, dict):
            continue
        source = item.get("source")
        if not isinstance(source, dict):
            continue
        command = source.get("download_command")
        url = source.get("remote_url")
        rel = item.get("path")
        if command:
            lines.append(str(command))
        elif url and rel:
            lines.append(f'curl -L "{url}" -o "{rel}"')
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    out_path.chmod(0o755)


def write_compare_to_s3_sh(record: dict[str, Any], out_path: Path) -> None:
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", ""]
    target = record.get("replay", {}).get("target_output", {})
    if isinstance(target, dict):
        published = target.get("published")
        if isinstance(published, dict):
            uri = published.get("remote_url") or published.get("s3_url") or published.get("uri")
            rel = target.get("path")
            if uri and rel:
                lines.append(f'aws s3 cp "{uri}" - | md5sum')
                lines.append(f'md5sum "{rel}"')
                lines.append("")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    out_path.chmod(0o755)
