from __future__ import annotations

import csv
import gzip
import io
import json
from pathlib import Path
import tarfile

from geneset_extractors.cli import main


def _write_gzip(path: Path, text: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write(text)


def _write_miniml(path: Path) -> None:
    sample_xml = []
    samples = [
        ("GSM1", "case_1", "COVID-19"),
        ("GSM2", "case_2", "COVID-19"),
        ("GSM3", "control_1", "Healthy"),
        ("GSM4", "control_2", "Healthy"),
        ("GSM5", "excluded_1", "Convalescent"),
    ]
    for accession, title, state in samples:
        sample_xml.append(
            f"""
  <Sample iid="{accession}">
    <Title>{title}</Title>
    <Accession database="GEO">{accession}</Accession>
    <Channel position="1">
      <Source>PBMC</Source>
      <Organism taxid="9606">Homo sapiens</Organism>
      <Characteristics tag="disease state">{state}</Characteristics>
      <Characteristics tag="cell type">PBMC</Characteristics>
    </Channel>
  </Sample>"""
        )
    payload = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<MINiML xmlns="http://www.ncbi.nlm.nih.gov/geo/info/MINiML" version="0.5.0">\n'
        + "\n".join(sample_xml)
        + "\n</MINiML>\n"
    )
    data = payload.encode("utf-8")
    info = tarfile.TarInfo("GSETEST_family.xml")
    info.size = len(data)
    with tarfile.open(path, "w:gz") as archive:
        archive.addfile(info, io.BytesIO(data))


def _make_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    counts = tmp_path / "counts.tsv.gz"
    miniml = tmp_path / "family.xml.tgz"
    annotation = tmp_path / "annotation.tsv.gz"
    _write_gzip(
        counts,
        "gene_id\tcase_1\tcase_2\tcontrol_1\tcontrol_2\texcluded_1\n"
        "ENSG1\t100\t110\t10\t12\t30\n"
        "ENSG2\t8\t10\t90\t100\t20\n"
        "ENSG3\t40\t42\t41\t39\t40\n",
    )
    _write_miniml(miniml)
    _write_gzip(
        annotation,
        "EnsemblGeneID\tSymbol\nENSG1\tUP1\nENSG2\tDOWN1\nENSG3\tFLAT1\n",
    )
    return counts, miniml, annotation


def test_geo_bulk_prepare_and_de_provenance_chain(tmp_path: Path) -> None:
    counts, miniml, annotation = _make_inputs(tmp_path)
    prepared = tmp_path / "prepared"
    rc = main(
        [
            "workflows",
            "geo_bulk_prepare",
            "--counts_file",
            str(counts),
            "--miniml_file",
            str(miniml),
            "--annotation_file",
            str(annotation),
            "--out_dir",
            str(prepared),
            "--study_id",
            "GSETEST",
            "--group_characteristic",
            "disease state",
            "--condition_a_values",
            "COVID-19",
            "--condition_b_values",
            "Healthy",
            "--counts_source_url",
            "https://example.org/GSETEST_counts.tsv.gz",
            "--miniml_source_url",
            "https://example.org/GSETEST_family.xml.tgz",
            "--annotation_source_url",
            "https://example.org/annotation.tsv.gz",
            "--landing_page_url",
            "https://example.org/GSETEST",
        ]
    )
    assert rc == 0

    with (prepared / "sample_metadata.tsv").open(encoding="utf-8", newline="") as handle:
        metadata_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["sample_id"] for row in metadata_rows] == ["case_1", "case_2", "control_1", "control_2"]
    assert {row["condition"] for row in metadata_rows} == {"case", "control"}

    counts_header = (prepared / "counts.tsv").read_text(encoding="utf-8").splitlines()[0]
    assert counts_header == "gene_id\tcase_1\tcase_2\tcontrol_1\tcontrol_2"

    de_out = tmp_path / "de"
    rc = main(
        [
            "workflows",
            "rna_de_prepare",
            "--modality",
            "bulk",
            "--counts_tsv",
            str(prepared / "counts.tsv"),
            "--matrix_orientation",
            "gene_by_sample",
            "--feature_id_column",
            "gene_id",
            "--sample_id_column",
            "sample_id",
            "--sample_metadata_tsv",
            str(prepared / "sample_metadata.tsv"),
            "--group_column",
            "condition",
            "--comparison_mode",
            "condition_a_vs_b",
            "--condition_a",
            "case",
            "--condition_b",
            "control",
            "--feature_mapping_tsv",
            str(prepared / "feature_mapping.tsv"),
            "--feature_mapping_from_column",
            "source_feature_id",
            "--feature_mapping_to_column",
            "gene_symbol",
            "--backend",
            "lightweight",
            "--out_dir",
            str(de_out),
            "--upstream_provenance_graph_json",
            str(prepared / "geo_bulk_inputs.provenance_graph.json"),
        ]
    )
    assert rc == 0

    graph = json.loads((de_out / "deg_long.provenance_graph.json").read_text(encoding="utf-8"))["deg_long"]
    methods = {node.get("analysis", {}).get("script_url", "") for node in graph["nodes"] if node.get("type") == "AnalysisType"}
    names = {node.get("name") for node in graph["nodes"] if node.get("type") == "AnalysisType"}
    urls = {
        value
        for node in graph["nodes"]
        if node.get("type") == "File"
        for value in (node.get("dcc_url"), node.get("drc_url"), node.get("c2m2_properties", {}).get("local_id"))
    }
    assert "prepare_GSETEST" in names
    assert "prepare_deg_long" in names
    assert "https://example.org/GSETEST_counts.tsv.gz" in urls
    assert methods
    node_id_by_name = {node.get("name"): node.get("id") for node in graph["nodes"]}
    prepare_id = node_id_by_name["prepare_GSETEST"]
    de_id = node_id_by_name["prepare_deg_long"]
    prepared_counts_id = node_id_by_name["counts.tsv"]
    edge_pairs = {(edge["source"], edge["target"]) for edge in graph["edges"]}
    assert (prepare_id, prepared_counts_id) in edge_pairs
    assert (prepared_counts_id, de_id) in edge_pairs


def test_geo_bulk_prepare_can_match_standardized_counts_by_gsm_accession(tmp_path: Path) -> None:
    _, miniml, annotation = _make_inputs(tmp_path)
    counts = tmp_path / "accession_counts.tsv.gz"
    _write_gzip(
        counts,
        "gene_id\tGSM1\tGSM2\tGSM3\tGSM4\tGSM5\n"
        "ENSG1\t100\t110\t10\t12\t30\n"
        "ENSG2\t8\t10\t90\t100\t20\n",
    )
    prepared = tmp_path / "prepared_accessions"
    rc = main(
        [
            "workflows",
            "geo_bulk_prepare",
            "--counts_file",
            str(counts),
            "--miniml_file",
            str(miniml),
            "--annotation_file",
            str(annotation),
            "--out_dir",
            str(prepared),
            "--study_id",
            "GSETEST",
            "--sample_id_field",
            "accession",
            "--group_characteristic",
            "disease state",
            "--condition_a_values",
            "COVID-19",
            "--condition_b_values",
            "Healthy",
        ]
    )
    assert rc == 0
    counts_header = (prepared / "counts.tsv").read_text(encoding="utf-8").splitlines()[0]
    assert counts_header == "gene_id\tGSM1\tGSM2\tGSM3\tGSM4"
    with (prepared / "sample_metadata.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["sample_id"] for row in rows] == ["GSM1", "GSM2", "GSM3", "GSM4"]
    assert [row["geo_title"] for row in rows] == ["case_1", "case_2", "control_1", "control_2"]


def test_geo_bulk_prepare_can_group_by_sample_title(tmp_path: Path) -> None:
    counts, miniml, annotation = _make_inputs(tmp_path)
    prepared = tmp_path / "prepared_title_groups"
    rc = main(
        [
            "workflows",
            "geo_bulk_prepare",
            "--counts_file",
            str(counts),
            "--miniml_file",
            str(miniml),
            "--annotation_file",
            str(annotation),
            "--out_dir",
            str(prepared),
            "--study_id",
            "GSETEST",
            "--group_characteristic",
            "__title__",
            "--condition_a_values",
            "case_1,case_2",
            "--condition_b_values",
            "control_1,control_2",
        ]
    )
    assert rc == 0
    with (prepared / "sample_metadata.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["condition"] for row in rows] == ["case", "case", "control", "control"]
    assert [row["geo_group_value"] for row in rows] == ["case_1", "case_2", "control_1", "control_2"]
