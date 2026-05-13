import json
from pathlib import Path

import pytest

from geneset_extractors.core.metadata import make_metadata
from geneset_extractors.core.provenance import REPO_URL, activate_runtime_context, build_analysis_node, build_file_node, mirror_graph_payload
from geneset_extractors.hashing import sha256_file
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


def test_build_file_node_can_mirror_local_identity_to_remote_prefix(tmp_path: Path):
    local_file = tmp_path / "outputs" / "geneset.meta.json"
    local_file.parent.mkdir(parents=True, exist_ok=True)
    local_file.write_text('{"ok": true}\n', encoding="utf-8")

    local_node = build_file_node(
        {"path": str(local_file), "local_path": str(local_file), "role": "metadata_json"},
        {},
    )
    mirrored_node = build_file_node(
        {"path": str(local_file), "local_path": str(local_file), "role": "metadata_json"},
        {},
        mirror_local_prefix=str(tmp_path),
        mirror_remote_prefix="s3://bucket/published",
    )

    assert local_node["c2m2_properties"]["md5"] == mirrored_node["c2m2_properties"]["md5"]
    assert mirrored_node["c2m2_properties"]["local_id"] == "s3://bucket/published/outputs/geneset.meta.json"
    assert mirrored_node["dcc_url"] == "s3://bucket/published/outputs/geneset.meta.json"
    assert mirrored_node["drc_url"] == "s3://bucket/published/outputs/geneset.meta.json"
    assert mirrored_node["id"] != local_node["id"]


def test_mirror_graph_payload_rewrites_upstream_file_nodes_and_commands(tmp_path: Path):
    local_file = tmp_path / "workflow" / "deg_long.tsv"
    local_file.parent.mkdir(parents=True, exist_ok=True)
    local_file.write_text("x\n", encoding="utf-8")
    local_uri = local_file.resolve().as_uri()
    payload = {
        "deg_long": {
            "nodes": [
                {
                    "id": "file:deg_tsv:old",
                    "type": "File",
                    "name": "deg_long.tsv",
                    "description": "File node for deg_long.tsv (role: deg_tsv)",
                    "dcc_url": local_uri,
                    "drc_url": local_uri,
                    "c2m2_properties": {
                        "filename": "deg_long.tsv",
                        "local_id": str(local_file.resolve()),
                        "size_in_bytes": 2,
                        "md5": "ndTkYSaMgDT1yFZOFVxnpg==",
                        "_uuid": "old",
                        "persistent_id": "old-persistent",
                    },
                },
                {
                    "id": "analysis:rna_de_prepare:old",
                    "type": "AnalysisType",
                    "name": "prepare_deg_long",
                    "description": "Analysis step",
                    "dcc_url": REPO_URL,
                    "drc_url": REPO_URL,
                    "c2m2_properties": {"id": "analysis:rna_de_prepare:old", "name": "prepare_deg_long", "description": "Analysis step", "synonyms": []},
                    "analysis": {
                        "script_url": REPO_URL,
                        "command": f"python script.py --counts_tsv {local_file.resolve()}",
                        "parameters": {"path": str(local_file.resolve())},
                        "environment": {"entrypoint": str(local_file.resolve())},
                    },
                },
            ],
            "edges": [
                {
                    "id": "file:deg_tsv:old_to_analysis:rna_de_prepare:old",
                    "source": "file:deg_tsv:old",
                    "target": "analysis:rna_de_prepare:old",
                    "label": "data input",
                    "description": "edge",
                }
            ],
        }
    }
    mirrored = mirror_graph_payload(
        payload,
        str(tmp_path),
        "s3://bucket/published",
    )
    graph = mirrored["deg_long"]
    file_node = next(node for node in graph["nodes"] if node["type"] == "File")
    analysis_node = next(node for node in graph["nodes"] if node["type"] == "AnalysisType")
    assert file_node["c2m2_properties"]["local_id"] == "s3://bucket/published/workflow/deg_long.tsv"
    assert file_node["dcc_url"] == "s3://bucket/published/workflow/deg_long.tsv"
    assert "s3://bucket/published/workflow/deg_long.tsv" in analysis_node["analysis"]["command"]
    assert analysis_node["analysis"]["parameters"]["path"] == "s3://bucket/published/workflow/deg_long.tsv"
    assert analysis_node["analysis"]["environment"]["entrypoint"] == "s3://bucket/published/workflow/deg_long.tsv"
    assert graph["edges"][0]["source"] == file_node["id"]


def test_build_analysis_node_normalizes_local_cli_paths_to_portable_forms():
    node = build_analysis_node(
        analysis_id="analysis:test",
        method="rna_deg_multi",
        name="generate_test",
        description="test",
        parameters={},
        command=[
            "/home/ryank/software/geneset_extractors/dig-gene-set-extractors/src/geneset_extractors/cli.py",
            "convert",
            "rna_deg_multi",
            "--deg_tsv",
            "s3://bucket/file.tsv",
        ],
        entrypoint="/home/ryank/software/geneset_extractors/dig-gene-set-extractors/src/geneset_extractors/cli.py convert rna_deg_multi",
        repo_url=REPO_URL,
    )
    assert node["analysis"]["command"].startswith("python -m geneset_extractors.cli convert rna_deg_multi ")
    assert node["analysis"]["environment"]["entrypoint"] == "geneset-extractors convert rna_deg_multi"


def test_mirror_graph_payload_preserves_file_identity_consistent_with_build_file_node(tmp_path: Path):
    local_file = tmp_path / "workflow" / "deg_long.tsv"
    local_file.parent.mkdir(parents=True, exist_ok=True)
    local_file.write_text("x\n", encoding="utf-8")
    record = {
        "path": str(local_file),
        "local_path": str(local_file),
        "role": "deg_tsv",
        "sha256": sha256_file(local_file),
        "size_bytes": local_file.stat().st_size,
    }
    original = build_file_node(
        record,
        {},
    )
    payload = {
        "deg_long": {
            "nodes": [original],
            "edges": [],
        }
    }
    mirrored_payload = mirror_graph_payload(payload, str(tmp_path), "s3://bucket/published")
    mirrored_node = mirrored_payload["deg_long"]["nodes"][0]
    rebuilt_node = build_file_node(
        record,
        {},
        mirror_local_prefix=str(tmp_path),
        mirror_remote_prefix="s3://bucket/published",
    )
    assert mirrored_node["id"] == rebuilt_node["id"]


def test_make_metadata_stores_effective_command_and_preserves_observed_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    input_file = tmp_path / "deg_long.tsv"
    input_file.write_text("x\n", encoding="utf-8")
    monkeypatch.setattr(
        __import__("sys"),
        "argv",
        [
            "/home/ryank/software/geneset_extractors/dig-gene-set-extractors/src/geneset_extractors/cli.py",
            "convert",
            "rna_deg_multi",
            "--score_mode",
            "auto",
            "--postprocess_mode",
            "harmonizome",
        ],
    )
    activate_runtime_context("rna_deg_multi")
    meta = make_metadata(
        converter_name="rna_deg_multi",
        parameters={
            "signature_name": "AB1",
            "comparison_label": "age30_20",
            "score_mode": "signed_neglog10padj",
            "score_mode_requested": "signed_neglog10padj",
            "padj_max": 0.05,
            "emit_small_gene_sets": True,
        },
        data_type="rna_seq",
        assay="bulk",
        organism="human",
        genome_build="hg38",
        files=[{"path": str(input_file), "local_path": str(input_file), "role": "deg_tsv", "sha256": "abc", "size_bytes": 2}],
        command_io={
            "deg_tsv": str(input_file),
            "out_dir": str(tmp_path / "out"),
            "organism": "human",
            "genome_build": "hg38",
        },
        gene_annotation={"mode": "none"},
        weights={"weight_type": "signed"},
        summary={"n_genes": 1},
    )
    execution = meta["converter"]["execution"]
    command = execution["command"]
    observed = execution["observed_command"]
    assert command[:4] == [__import__("sys").executable, "-m", "geneset_extractors.cli", "convert"]
    assert "--deg_tsv" in command and str(input_file) in command
    assert "--out_dir" in command and str(tmp_path / "out") in command
    assert "--organism" in command and "human" in command
    assert "--genome_build" in command and "hg38" in command
    assert "--score_mode" in command and "signed_neglog10padj" in command
    assert "--padj_max" in command and "0.05" in command
    assert observed[0].endswith("cli.py")
    assert "--score_mode" in observed and "auto" in observed
