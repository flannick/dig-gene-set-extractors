#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Offline installer for omics2geneset.

Usage:
  scripts/install_offline.sh --wheelhouse /path/to/wheels [--python /path/to/python] [--no-dev] [--wheel-install]

Options:
  --wheelhouse PATH   Directory containing pre-downloaded wheels (required).
  --python PATH       Python interpreter to use (default: python).
  --no-dev            Install package without [dev] extras.
  --wheel-install     Install non-editable wheel-style package (default is editable install).

Notes:
  - This script bootstraps build prerequisites (setuptools, wheel, build) first.
  - It uses --no-index and --no-build-isolation for air-gapped environments.
USAGE
}

WHEELHOUSE=""
PY_BIN="python"
WITH_DEV=1
EDITABLE=1
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --wheelhouse)
      WHEELHOUSE="${2:-}"
      shift 2
      ;;
    --python)
      PY_BIN="${2:-}"
      shift 2
      ;;
    --no-dev)
      WITH_DEV=0
      shift
      ;;
    --wheel-install)
      EDITABLE=0
      shift
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

if [[ -z "$WHEELHOUSE" ]]; then
  echo "--wheelhouse is required" >&2
  usage >&2
  exit 2
fi

if [[ ! -d "$WHEELHOUSE" ]]; then
  echo "Wheelhouse not found: $WHEELHOUSE" >&2
  exit 2
fi

if [[ "$WITH_DEV" -eq 1 ]]; then
  "$PY_BIN" "$SCRIPT_DIR/check_wheelhouse.py" --wheelhouse "$WHEELHOUSE" --with-dev
else
  "$PY_BIN" "$SCRIPT_DIR/check_wheelhouse.py" --wheelhouse "$WHEELHOUSE"
fi

TARGET="."
if [[ "$WITH_DEV" -eq 1 ]]; then
  TARGET=".[dev]"
fi

echo "Bootstrapping build prerequisites from wheelhouse..."
"$PY_BIN" -m pip install --no-index --find-links "$WHEELHOUSE" setuptools wheel build
"$PY_BIN" - <<'PY'
import importlib.util
import sys
if importlib.util.find_spec("wheel") is None:
    print("wheel package is still missing after bootstrap; cannot continue with --no-build-isolation install", file=sys.stderr)
    raise SystemExit(2)
PY

echo "Installing omics2geneset from local source..."
if [[ "$EDITABLE" -eq 1 ]]; then
  "$PY_BIN" -m pip install --no-index --find-links "$WHEELHOUSE" -e "$TARGET" --no-build-isolation
else
  "$PY_BIN" -m pip install --no-index --find-links "$WHEELHOUSE" "$TARGET" --no-build-isolation
fi

echo "Install complete. Verifying CLI..."
"$PY_BIN" -m omics2geneset.cli list >/dev/null
echo "Verification passed."
