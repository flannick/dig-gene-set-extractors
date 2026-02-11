import json
from pathlib import Path

import pytest

from omics2geneset.core.validate import validate_metadata_schema


def test_metadata_missing_required_fails(tmp_path: Path):
    meta = tmp_path / "geneset.meta.json"
    schema = Path("src/omics2geneset/schemas/geneset_metadata.schema.json")
    meta.write_text(json.dumps({"schema_version": "1.0.0"}), encoding="utf-8")
    with pytest.raises(Exception):
        validate_metadata_schema(meta, schema)
