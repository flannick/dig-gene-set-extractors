import json
from pathlib import Path

import pytest

from geneset_extractors.core.provenance import build_file_node
from geneset_extractors.core.validate import validate_metadata_schema, validate_provenance_schema


def test_metadata_missing_required_fails(tmp_path: Path):
    meta = tmp_path / "geneset.meta.json"
    schema = Path("src/geneset_extractors/schemas/geneset_metadata.schema.json")
    meta.write_text(json.dumps({"schema_version": "1.0.0"}), encoding="utf-8")
    with pytest.raises(Exception):
        validate_metadata_schema(meta, schema)


def test_provenance_missing_required_fails(tmp_path: Path):
    provenance = tmp_path / "geneset.provenance.json"
    schema = Path("src/geneset_extractors/schemas/geneset_provenance.schema.json")
    provenance.write_text(json.dumps({"graph1": {"nodes": [], "edges": []}}), encoding="utf-8")
    with pytest.raises(Exception):
        validate_provenance_schema(provenance, schema)


def test_build_file_node_raises_when_md5_cannot_be_computed(tmp_path: Path):
    missing = tmp_path / "missing.tsv"
    with pytest.raises(ValueError, match="Unable to compute MD5"):
        build_file_node({"path": str(missing), "local_path": str(missing), "role": "deg_tsv"}, {})
