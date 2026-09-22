"""Verified import of externally generated GMTs without scientific reprocessing."""
from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any

from geneset_extractors.core.metadata import input_file_record, make_metadata, write_metadata
from geneset_extractors.core.validate import validate_output_dir
from geneset_extractors.hashing import sha256_file


def _gmt_summary(path: Path) -> tuple[int, int]:
    names: set[str] = set()
    genes: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or not fields[0].strip():
                raise ValueError(f"malformed GMT row {line_number}: expected name, description, and at least one gene")
            if fields[0] in names:
                raise ValueError(f"duplicate GMT set name at row {line_number}: {fields[0]}")
            names.add(fields[0])
            genes.update(gene for gene in fields[2:] if gene)
    if not names:
        raise ValueError("GMT contains no gene sets")
    return len(names), len(genes)


def import_external_gmt(*, gmt: Path, out_dir: Path, source_record: Path, library_id: str, model_id: str, display_name: str, description: str, expected_sha256: str | None = None) -> dict[str, Any]:
    """Copy one external GMT unchanged and emit truthful metadata/provenance."""
    gmt = gmt.resolve()
    if not gmt.is_file():
        raise ValueError(f"external GMT does not exist: {gmt}")
    source = json.loads(source_record.read_text(encoding="utf-8"))
    if not isinstance(source, dict):
        raise ValueError("source record must be a JSON object")
    required = ("name", "uri_or_identifier", "release", "license", "access_restrictions", "organism", "genome_build", "assay", "data_type")
    missing = [key for key in required if not str(source.get(key, "")).strip()]
    if missing:
        raise ValueError(f"source record missing required fields: {', '.join(missing)}")
    set_count, gene_count = _gmt_summary(gmt)
    out_dir.mkdir(parents=True, exist_ok=True)
    destination = out_dir / "genesets.gmt"
    shutil.copyfile(gmt, destination)
    source_checksum = sha256_file(gmt)
    if expected_sha256 and source_checksum != expected_sha256.removeprefix("sha256:"):
        raise ValueError(f"external GMT checksum mismatch: expected {expected_sha256}, observed sha256:{source_checksum}")
    if sha256_file(destination) != source_checksum:
        raise RuntimeError("copied GMT checksum differs from source GMT")
    source_file = input_file_record(
        gmt,
        "external_precomputed_gmt",
        canonical_uri=str(source["uri_or_identifier"]),
        version=str(source["release"]),
        license=str(source["license"]),
        access_level=str(source["access_restrictions"]),
        provider=str(source["name"]),
    )
    metadata = make_metadata(
        converter_name="external_precomputed_import",
        parameters={
            "library_id": library_id,
            "model_id": model_id,
            "source_release": str(source["release"]),
            "generation_status": "external_precomputed_incomplete_code",
            "source_gmt_sha256": source_checksum,
        },
        data_type=str(source["data_type"]),
        assay=str(source["assay"]),
        organism=str(source["organism"]),
        genome_build=str(source["genome_build"]),
        files=[source_file],
        gene_annotation={"mode": "none", "source": "external precomputed GMT", "gene_id_field": "gene_symbol"},
        weights={"weight_type": "unweighted", "normalization": {}, "aggregation": "external_precomputed"},
        # The external artifact is already a set of gene identifiers: every
        # observed identifier is retained, with no repository-side mapping or
        # selection. These fields satisfy the standard metadata schema without
        # implying that DIG performed the original scientific analysis.
        summary={
            "n_input_features": gene_count,
            "n_genes": gene_count,
            "n_features_assigned": gene_count,
            "fraction_features_assigned": 1.0,
            "n_gene_sets": set_count,
            "generation_status": "external_precomputed_incomplete_code",
        },
        output_files=[{"path": "genesets.gmt", "role": "gmt"}],
        gmt={"written": True, "path": "genesets.gmt", "prefer_symbol": True, "min_genes": 0, "max_genes": 0, "plans": []},
        gene_set_description=description,
        command_io={"input": str(gmt), "out_dir": str(out_dir)},
    )
    metadata["gene_set"]["name"] = display_name  # type: ignore[index]
    metadata["gene_set"]["primary_artifact"] = {"path": "genesets.gmt", "role": "gmt"}  # type: ignore[index]
    # ``make_metadata`` supplies the normal extracted-program TSV by default.
    # An external import has no such derived artifact: its verified GMT is the
    # only data output, so retaining that default would make provenance hash a
    # nonexistent geneset.tsv.
    metadata["output"]["files"] = [  # type: ignore[index]
        item for item in metadata["output"]["files"]  # type: ignore[index]
        if item.get("path") != "geneset.tsv"
    ]
    metadata["converter"]["code"].update({  # type: ignore[index]
        "repo_url": str(source["uri_or_identifier"]),
        "module": None,
        "script_url": source.get("documentation"),
        "description": "Externally generated gene sets imported unchanged after checksum verification; complete regeneration code was not supplied.",
    })
    metadata["converter"]["execution"].update({  # type: ignore[index]
        "mode": "external_precomputed_import",
        "entrypoint": "geneset-extractors external-import",
        "notes": "Repository activity validates and imports the supplied GMT; it does not reproduce external scientific generation.",
    })
    metadata["external_import"] = {"source_name": source["name"], "regeneration_status": "incomplete_code", "source_gmt_sha256": source_checksum}
    write_metadata(out_dir / "geneset.meta.json", metadata)
    validate_output_dir(out_dir, Path(__file__).parent / "schemas" / "geneset_metadata.schema.json")
    return {"out_dir": str(out_dir), "model_id": model_id, "sha256": source_checksum, "n_gene_sets": set_count, "n_genes": gene_count}
