import json
from pathlib import Path

import yaml

from geneset_extractors.core.dapper_provenance import (
    DAPPER_RELEASE,
    DAPPER_SCHEMA_REVISION,
    _compute_id,
    build_dapper_provenance,
)
from geneset_extractors.core.metadata import (
    DAPPER_PROVENANCE_FILENAME,
    LEGACY_PROVENANCE_FILENAME,
    input_file_record,
    make_metadata,
    write_metadata,
)
from geneset_extractors.core.validate import validate_output_dir


def test_metadata_write_emits_renamed_legacy_and_dapper_provenance(tmp_path: Path):
    source = tmp_path / "source.tsv"
    source.write_text("gene_id\tscore\nGENE1\t1\n", encoding="utf-8")
    geneset = tmp_path / "geneset.tsv"
    geneset.write_text("gene_id\tscore\nGENE1\t1\n", encoding="utf-8")
    metadata = make_metadata(
        converter_name="toy_converter",
        parameters={"term_prefix": "Toy"},
        data_type="expression",
        assay="bulk_rna",
        organism="human",
        genome_build="hg38",
        files=[input_file_record(source, "source_tsv")],
        gene_annotation={"mode": "none", "source": "toy", "gene_id_field": "symbol"},
        weights={"weight_type": "score", "normalization": {}, "aggregation": "none"},
        summary={
            "n_input_features": 1,
            "n_genes": 1,
            "n_features_assigned": 1,
            "fraction_features_assigned": 1.0,
            "n_sets_emitted": 1,
        },
    )

    write_metadata(tmp_path / "geneset.meta.json", metadata)

    legacy_path = tmp_path / LEGACY_PROVENANCE_FILENAME
    dapper_path = tmp_path / DAPPER_PROVENANCE_FILENAME
    assert legacy_path.exists()
    assert dapper_path.exists()
    assert not (tmp_path / "geneset.provenance.json").exists()

    meta_payload = json.loads((tmp_path / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta_payload["provenance"] == {
        "path": LEGACY_PROVENANCE_FILENAME,
        "dapper_path": DAPPER_PROVENANCE_FILENAME,
        "focus_node_id": meta_payload["gene_set"]["id"],
    }
    dapper_payload = yaml.safe_load(dapper_path.read_text(encoding="utf-8"))
    assert {"c2m2_files", "activities", "gene_sets", "used_edges", "was_generated_by_edges"}.issubset(dapper_payload)
    all_node_ids = {
        node["id"]
        for bucket in ("c2m2_files", "activities", "gene_sets")
        for node in dapper_payload[bucket]
    }
    assert all(node_id.startswith("dapper:") for node_id in all_node_ids)
    assert dapper_payload["used_edges"][0]["predicate"] == "prov:used"
    assert dapper_payload["was_generated_by_edges"][0]["predicate"] == "prov:wasGeneratedBy"
    assert dapper_payload["used_edges"][0]["subject"] in all_node_ids
    assert dapper_payload["used_edges"][0]["object"] in all_node_ids
    assert dapper_payload["was_generated_by_edges"][0]["subject"] in all_node_ids
    assert dapper_payload["was_generated_by_edges"][0]["object"] in all_node_ids
    assert validate_output_dir(
        tmp_path,
        Path("src/geneset_extractors/schemas/geneset_metadata.schema.json"),
    ) == {"mode": "single", "n_groups": 1}


def test_dapper_0_2_routes_non_c2m2_file_nodes_to_generic_file_bucket():
    """Dapper 0.2 distinguishes generic files from C2M2-registered files."""
    legacy_payload = {
        "toy": {
            "nodes": [
                {
                    "id": "source-file",
                    "type": "File",
                    "name": "source.tsv",
                    "filename": "source.tsv",
                    "location": "https://example.test/source.tsv",
                    "sha256": "a" * 64,
                    "mime_type": "text/tab-separated-values",
                },
                {"id": "operation", "type": "AnalysisType", "name": "extract"},
                {"id": "result", "type": "GeneSet", "name": "toy"},
            ],
            "edges": [
                {"source": "source-file", "target": "operation", "label": "data input"},
                {"source": "operation", "target": "result", "label": "data output"},
            ],
        }
    }
    payload = build_dapper_provenance(legacy_payload, {})

    assert DAPPER_RELEASE == "0.2.0-a0"
    assert DAPPER_SCHEMA_REVISION == "af9f391fdcc64a0d1bc3a4f3073c0fff6a55e968"
    assert "files" in payload
    assert "c2m2_files" not in payload
    file_node = payload["files"][0]
    assert file_node["id"].startswith("dapper:File.")
    assert file_node["location"] == "https://example.test/source.tsv"
    node_ids = {
        node["id"]
        for bucket in ("files", "activities", "gene_sets")
        for node in payload[bucket]
    }
    assert payload["used_edges"][0]["object"] in node_ids
    assert payload["was_generated_by_edges"][0]["subject"] in node_ids


def test_dapper_0_2_file_identity_ignores_location_and_preserves_literal_scalars():
    """Dapper 0.2 does not rewrite a literal scalar that resembles an ID."""
    first = {
        "id": "external-file-id",
        "name": "external-file-id",
        "filename": "source.tsv",
        "location": "/one/work/source.tsv",
    }
    second = dict(first, location="s3://example/source.tsv")

    first_id = _compute_id(first, "File", first["id"])
    # Simulate a second identity check after the node has been minted. The
    # scalar name remains external metadata, not a DAPPER node reference.
    first["id"] = first_id
    assert _compute_id(first, "File", first["id"]) == first_id
    assert _compute_id(second, "File", second["id"]) == first_id
