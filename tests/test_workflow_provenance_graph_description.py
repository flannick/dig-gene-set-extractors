import json
from pathlib import Path

from geneset_extractors.workflows.gtex_runtime_common import write_workflow_provenance_graph


def _analysis_desc(graph_path: Path) -> str:
    payload = json.loads(graph_path.read_text(encoding="utf-8"))
    graph = list(payload.values())[0]
    analysis = [n for n in graph["nodes"] if n["type"] == "AnalysisType"]
    assert len(analysis) == 1
    return analysis[0]["description"]


def test_custom_description_flows_into_operation_node(tmp_path: Path):
    inp = tmp_path / "raw.tsv"
    inp.write_text("a\tb\n1\t2\n", encoding="utf-8")
    out = tmp_path / "prepared.tsv"
    out.write_text("x\ty\n3\t4\n", encoding="utf-8")
    graph_path = write_workflow_provenance_graph(
        workflow_name="ptm_prepare_public",
        module_name="geneset_extractors.workflows.ptm_prepare_public",
        output_dir=tmp_path,
        focus_output_path=out,
        output_paths=[(out, "prepared")],
        input_paths=[(inp, "raw")],
        parameters={"ptm_type": "phospho"},
        description="CUSTOM CPTAC DESCRIPTION",
    )
    assert _analysis_desc(graph_path) == "CUSTOM CPTAC DESCRIPTION"


def test_default_description_unchanged(tmp_path: Path):
    inp = tmp_path / "raw.tsv"
    inp.write_text("a\n1\n", encoding="utf-8")
    out = tmp_path / "prepared.tsv"
    out.write_text("x\n3\n", encoding="utf-8")
    graph_path = write_workflow_provenance_graph(
        workflow_name="gtex_age_binned",
        module_name="geneset_extractors.workflows.gtex_age_binned",
        output_dir=tmp_path,
        focus_output_path=out,
        output_paths=[(out, "prepared")],
        input_paths=[(inp, "raw")],
        parameters={},
    )
    assert "GTEx differential expression results" in _analysis_desc(graph_path)
