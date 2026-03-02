"""CLI compatibility entry point for the neutral geneset_extractors namespace."""

from __future__ import annotations

from omics2geneset.cli import build_parser, main

__all__ = ["build_parser", "main"]


if __name__ == "__main__":
    raise SystemExit(main())

