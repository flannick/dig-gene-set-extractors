import hashlib
import json
from pathlib import Path
import subprocess
import sys

import yaml

from geneset_extractors.core.dapper_provenance import validate_dapper_document
from geneset_extractors.core.provenance_convert import deduplicate_legacy_provenance


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "geneset_extractors.cli", *args],
        capture_output=True,
        text=True,
        env={**__import__("os").environ, "PYTHONPATH": "src"},
    )


def _write_old_dig_pair(directory: Path, name: str = "geneset.provenance.json") -> tuple[Path, Path]:
    directory.mkdir(parents=True, exist_ok=True)
    # Original DIG used this unwrapped {nodes, edges} layout before graph maps.
    legacy = {
        "nodes": [
            {"id": "source", "type": "File", "name": "source.tsv", "description": "source", "location": "source.tsv"},
            {"id": "extract", "type": "AnalysisType", "name": "extract", "description": "extract"},
            {"id": "result", "type": "GeneSet", "name": "old result", "description": "result"},
        ],
        "edges": [
            {"id": "input", "source": "source", "target": "extract", "label": "data input"},
            {"id": "output", "source": "extract", "target": "result", "label": "data output"},
        ],
    }
    metadata = {
        "gene_set": {"organism": "human", "assay": "bulk_rna", "data_type": "expression", "genome_build": "hg38"},
        "converter": {"parameters": {}}, "summary": {"n_genes": 1},
        "input": {"files": []}, "output": {"files": []},
    }
    provenance = directory / name
    meta = directory / "geneset.meta.json"
    provenance.write_text(json.dumps(legacy, indent=2) + "\n", encoding="utf-8")
    meta.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return provenance, meta


def test_convert_old_single_graph_discovers_metadata_and_preserves_inputs(tmp_path: Path):
    provenance, metadata = _write_old_dig_pair(tmp_path)
    before = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in (provenance, metadata)}
    result = _run("provenance", "convert", str(provenance))
    assert result.returncode == 0, result.stderr
    output = tmp_path / "geneset.provenance.dapper.yaml"
    assert "converted" in result.stdout and "summary converted=1 skipped=0 failed=0" in result.stdout
    assert output.exists()
    assert before == {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in (provenance, metadata)}
    document = yaml.safe_load(output.read_text(encoding="utf-8"))
    validate_dapper_document(document)
    assert document["used_edges"][0]["predicate"] == "prov:used"


def test_convert_directory_prefers_legacy_filename_and_recurses(tmp_path: Path):
    first, _ = _write_old_dig_pair(tmp_path / "one", "geneset.provenance.json")
    preferred, _ = _write_old_dig_pair(tmp_path / "two", "geneset.provenance.legacy.json")
    duplicate, _ = _write_old_dig_pair(tmp_path / "two", "geneset.provenance.json")
    third, _ = _write_old_dig_pair(tmp_path / "nested" / "three")
    discovered = _run("provenance", "discover", str(tmp_path), "--recursive")
    assert discovered.returncode == 0, discovered.stderr
    assert discovered.stdout.splitlines() == [str(third.resolve()), str(first.resolve()), str(preferred.resolve())]
    result = _run("provenance", "convert", str(tmp_path), "--recursive")
    assert result.returncode == 0, result.stderr
    assert "summary converted=3 skipped=0 failed=0" in result.stdout
    assert str(preferred) in result.stdout
    assert str(duplicate) not in result.stdout
    for provenance in (first, preferred, third):
        assert provenance.parent.joinpath("geneset.provenance.dapper.yaml").exists()


def test_convert_skips_existing_unless_overwritten_and_is_deterministic(tmp_path: Path):
    provenance, _ = _write_old_dig_pair(tmp_path, "geneset.provenance.legacy.json")
    output = tmp_path / "custom.dapper.yaml"
    first = _run("provenance", "convert", str(provenance), "--out", str(output))
    assert first.returncode == 0
    original = output.read_bytes()
    skipped = _run("provenance", "convert", str(provenance), "--out", str(output))
    assert skipped.returncode == 0 and "summary converted=0 skipped=1 failed=0" in skipped.stdout
    overwritten = _run("provenance", "convert", str(provenance), "--out", str(output), "--overwrite")
    assert overwritten.returncode == 0
    assert output.read_bytes() == original


def test_convert_reports_malformed_provenance_and_missing_or_invalid_metadata(tmp_path: Path):
    bad = tmp_path / "bad"; bad.mkdir()
    (bad / "geneset.provenance.json").write_text("{broken", encoding="utf-8")
    missing = _run("provenance", "convert", str(bad / "geneset.provenance.json"))
    assert missing.returncode == 1 and "metadata file is missing" in missing.stderr
    (bad / "geneset.meta.json").write_text("{}", encoding="utf-8")
    malformed = _run("provenance", "convert", str(bad / "geneset.provenance.json"))
    assert malformed.returncode == 1 and "not valid JSON" in malformed.stderr
    provenance, metadata = _write_old_dig_pair(tmp_path / "invalid")
    metadata.write_text("{}", encoding="utf-8")
    invalid = _run("provenance", "convert", str(provenance))
    assert invalid.returncode == 1 and "incompatible" in invalid.stderr


def test_recursive_conversion_continues_after_failure(tmp_path: Path):
    _write_old_dig_pair(tmp_path / "good")
    broken = tmp_path / "broken"; broken.mkdir()
    (broken / "geneset.provenance.legacy.json").write_text("{}", encoding="utf-8")
    (broken / "geneset.meta.json").write_text("{}", encoding="utf-8")
    result = _run("provenance", "convert", str(tmp_path), "--recursive")
    assert result.returncode == 1
    assert "summary converted=1 skipped=0 failed=1" in result.stdout
    assert (tmp_path / "good" / "geneset.provenance.dapper.yaml").exists()


def test_deduplicate_backs_up_and_removes_dapper_equivalent_file_nodes(tmp_path: Path):
    provenance, _ = _write_old_dig_pair(tmp_path)
    payload = json.loads(provenance.read_text(encoding="utf-8"))
    duplicate = dict(payload["nodes"][0])
    duplicate["id"] = "source-duplicate"
    payload["nodes"].append(duplicate)
    payload["edges"].append({"id": "input-duplicate", "source": "source-duplicate", "target": "extract", "label": "data input"})
    provenance.write_text(json.dumps(payload), encoding="utf-8")
    original = provenance.read_bytes()
    nodes, edges = deduplicate_legacy_provenance(provenance, overwrite=False)
    assert (nodes, edges) == (1, 1)
    assert provenance.with_name("geneset.provenance.duplicates.json").read_bytes() == original
    corrected = json.loads(provenance.read_text(encoding="utf-8"))
    assert len(corrected["nodes"]) == 3
    assert len(corrected["edges"]) == 2
