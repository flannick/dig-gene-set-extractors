from __future__ import annotations

import argparse
import pathlib
import re
import sys


def normalize_name(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def available_distributions(wheelhouse: pathlib.Path) -> set[str]:
    out: set[str] = set()
    for wheel in wheelhouse.glob("*.whl"):
        name = wheel.name.split("-")[0]
        out.add(normalize_name(name))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--wheelhouse", required=True)
    ap.add_argument("--with-dev", action="store_true")
    args = ap.parse_args()

    wheelhouse = pathlib.Path(args.wheelhouse)
    if not wheelhouse.is_dir():
        print(f"wheelhouse not found: {wheelhouse}", file=sys.stderr)
        return 2

    required = {"setuptools", "wheel", "build", "jsonschema"}
    if args.with_dev:
        required.add("pytest")

    available = available_distributions(wheelhouse)
    missing = sorted([pkg for pkg in required if normalize_name(pkg) not in available])

    if missing:
        print("wheelhouse missing required direct wheels:", file=sys.stderr)
        for pkg in missing:
            print(f"  - {pkg}", file=sys.stderr)
        return 2

    print("wheelhouse basic check passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
