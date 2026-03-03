"""Legacy CLI wrapper for backward compatibility."""

from __future__ import annotations

from geneset_extractors.cli import build_parser, main

__all__ = ["build_parser", "main"]


if __name__ == "__main__":
    raise SystemExit(main())
