from __future__ import annotations

import base64
import hashlib
from pathlib import Path


_CHUNK = 1024 * 1024


def _digest(path: Path, algorithm: str) -> bytes:
    h = hashlib.new(algorithm)
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(_CHUNK), b""):
            h.update(chunk)
    return h.digest()


def _digest_hex_stream(path: Path, algorithm: str) -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(_CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def md5_base64(path: Path) -> str:
    return base64.b64encode(_digest(path, "md5")).decode("ascii")


def sha256_hex(path: Path) -> str:
    if path.stat().st_size < _CHUNK:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    return _digest_hex_stream(path, "sha256")


def file_size(path: Path) -> int:
    return path.stat().st_size


def checksum_list(path: Path, existing_md5_base64: str | None = None) -> list[dict[str, str]]:
    return [
        {"algorithm": "md5", "encoding": "base64", "value": existing_md5_base64 or md5_base64(path)},
        {"algorithm": "sha256", "encoding": "hex", "value": sha256_hex(path)},
    ]

