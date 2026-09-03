import json
from pathlib import Path

import yaml

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
