from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path

from omics2geneset import __version__
from omics2geneset.hashing import sha256_file, stable_hash_object


def build_geneset_id(converter_name: str, file_hashes: list[str], params: dict[str, object]) -> str:
    return stable_hash_object(
        {
            "converter": converter_name,
            "file_hashes": sorted(file_hashes),
            "params_hash": stable_hash_object(params),
        }
    )[:24]


def write_metadata(path: str | Path, payload: dict[str, object]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def make_metadata(
    converter_name: str,
    parameters: dict[str, object],
    data_type: str,
    assay: str,
    organism: str,
    genome_build: str,
    files: list[dict[str, str]],
    gene_annotation: dict[str, object],
    weights: dict[str, object],
    summary: dict[str, object],
    program_extraction: dict[str, object] | None = None,
    output_files: list[dict[str, str]] | None = None,
    gmt: dict[str, object] | None = None,
) -> dict[str, object]:
    file_hashes = [f["sha256"] for f in files]
    geneset_id = build_geneset_id(converter_name, file_hashes, parameters)
    payload: dict[str, object] = {
        "schema_version": "1.0.0",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "geneset_id": geneset_id,
        "converter": {
            "name": converter_name,
            "version": __version__,
            "parameters": parameters,
            "code": {"git_commit": None},
        },
        "input": {
            "data_type": data_type,
            "assay": assay,
            "organism": organism,
            "genome_build": genome_build,
            "files": files,
        },
        "gene_annotation": gene_annotation,
        "weights": weights,
        "summary": summary,
    }
    if program_extraction is not None:
        payload["program_extraction"] = program_extraction
    if output_files is not None:
        payload["output"] = {"files": output_files}
    if gmt is not None:
        payload["gmt"] = gmt
    return payload


def input_file_record(path: str | Path, role: str) -> dict[str, str]:
    p = Path(path)
    return {"path": str(p), "role": role, "sha256": sha256_file(p)}
