#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Bootstrap build tools required for --no-build-isolation install paths.

Usage:
  scripts/bootstrap_build_tools.sh [--python /path/to/python] [--wheelhouse /path/to/wheels]

Behavior:
  - Installs/updates pip, setuptools, and wheel in the active environment.
  - If --wheelhouse is provided, installs from local wheels with --no-index.
USAGE
}

PY_BIN="python"
WHEELHOUSE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --python)
      PY_BIN="${2:-}"
      shift 2
      ;;
    --wheelhouse)
      WHEELHOUSE="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -n "$WHEELHOUSE" ]]; then
  if [[ ! -d "$WHEELHOUSE" ]]; then
    echo "Wheelhouse not found: $WHEELHOUSE" >&2
    exit 2
  fi
  "$PY_BIN" -m pip install --no-index --find-links "$WHEELHOUSE" -U pip setuptools wheel
else
  "$PY_BIN" -m pip install -U pip setuptools wheel
fi

echo "Build tools ready."
