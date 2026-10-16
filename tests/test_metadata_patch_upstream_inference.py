import json
from pathlib import Path

from geneset_extractors.core.metadata_patch import _infer_upstream_graph_from_metadata_payload


def _payload(ptm_matrix_path: str) -> dict:
    return {
        "converter": {"name": "ptm_site_matrix"},
        "input": {"files": [{"role": "ptm_matrix_tsv", "local_path": ptm_matrix_path}]},
    }


def test_infers_ptm_site_matrix_graph_when_sibling_exists(tmp_path: Path):
    (tmp_path / "ptm_matrix.tsv").write_text("x\n", encoding="utf-8")
    graph = tmp_path / "ptm_matrix.provenance_graph.json"
    graph.write_text("{}", encoding="utf-8")
    resolved = _infer_upstream_graph_from_metadata_payload(_payload(str(tmp_path / "ptm_matrix.tsv")))
    assert resolved == str(graph)


def test_returns_none_when_sibling_graph_absent(tmp_path: Path):
    (tmp_path / "ptm_matrix.tsv").write_text("x\n", encoding="utf-8")
    resolved = _infer_upstream_graph_from_metadata_payload(_payload(str(tmp_path / "ptm_matrix.tsv")))
    assert resolved is None


def test_other_converters_still_return_none(tmp_path: Path):
    payload = {"converter": {"name": "some_unrelated_converter"}, "input": {"files": []}}
    assert _infer_upstream_graph_from_metadata_payload(payload) is None
