from __future__ import annotations

import importlib.metadata as im
import platform
import subprocess
import sys
from pathlib import Path


def _run(cmd: list[str], cwd: Path | None = None) -> str | None:
    try:
        return subprocess.check_output(cmd, cwd=cwd, text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return None


def git_info(repo_path: Path) -> dict[str, object]:
    commit = _run(["git", "rev-parse", "HEAD"], cwd=repo_path)
    ref = _run(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=repo_path)
    status = _run(["git", "status", "--porcelain"], cwd=repo_path)
    return {
        "git_commit": commit or "unknown",
        "git_ref": ref,
        "git_dirty": bool(status),
    }


def _dist_version(name: str) -> str | None:
    try:
        return im.version(name)
    except im.PackageNotFoundError:
        return None


def python_environment() -> dict[str, object]:
    return {
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "package_version": _dist_version("geneset-extractors"),
    }


def pip_freeze_text() -> str | None:
    return _run([sys.executable, "-m", "pip", "freeze"])


def r_session_info_text() -> str | None:
    return _run(["Rscript", "-e", "writeLines(capture.output(sessionInfo()))"])


def r_package_version(pkg: str) -> str | None:
    return _run(["Rscript", "-e", f"cat(as.character(packageVersion('{pkg}')))"])

