from __future__ import annotations

from pathlib import Path
from typing import Any

def compare_record_to_mirrors(record: dict[str, Any], workdir: Path) -> dict[str, Any]:
    results: list[dict[str, str]] = []
    target = record.get("replay", {}).get("target_output", {})
    if isinstance(target, dict):
        published = target.get("published")
        if isinstance(published, dict):
            rel = str(target.get("path", "")).strip()
            mirror = str(published.get("remote_url") or published.get("s3_url") or published.get("uri") or "").strip()
            status = "missing_local"
            if rel and (workdir / rel).exists():
                status = "local_present"
            results.append({"id": str(target.get("id", "")), "mirror": mirror, "status": status})
    return {"n_targets": len(results), "results": results}
