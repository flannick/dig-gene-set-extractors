from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace

from geneset_extractors.cli import main
from geneset_extractors.workflows.igvf_perturbseq import run as run_igvf_perturbseq


GENES = [
    "TP53", "MYC", "EGFR", "KRAS", "BRCA1", "BRCA2", "PTEN", "RB1",
    "AKT1", "MTOR", "JAK2", "STAT3", "CDK4", "CDK6", "CCND1", "MDM2",
    "NF1", "APC", "VHL", "ATM", "CHEK2", "PIK3CA", "NRAS", "ERBB2",
]
PERTS = ["GENE_A_KO", "GENE_B_KO", "GENE_C_KO", "GENE_D_KO"]


def _write_perturbation_by_gene_matrix(path: Path) -> None:
    """Each perturbation spikes 3 genes up (+10) and 3 down (-10); the rest are 0."""
    rows: list[list[object]] = [["perturbation", *GENES]]
    for i, pert in enumerate(PERTS):
        values = {g: 0.0 for g in GENES}
        for g in GENES[i * 6 : i * 6 + 3]:
            values[g] = 10.0
        for g in GENES[i * 6 + 3 : i * 6 + 6]:
            values[g] = -10.0
        rows.append([pert, *[values[g] for g in GENES]])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        csv.writer(fh, delimiter="\t").writerows(rows)


def test_igvf_perturbseq_emits_signed_term_gene_and_provenance(tmp_path: Path):
    expr = tmp_path / "igvf_signatures.tsv"
    _write_perturbation_by_gene_matrix(expr)
    out_dir = tmp_path / "workflow"
    args = SimpleNamespace(
        expression_tsv=str(expr),
        mapping_file=None,
        out_dir=str(out_dir),
        organism="human",
        genome_build="hg38",
        orientation="perturbation_by_gene",
        gmt_name="gene_set_library_crisp.gmt",
        z_threshold=1.0,
        min_gmt_size=2,
        provenance_overlay_json=None,
        provenance_mirror_local_prefix=None,
        provenance_mirror_remote_prefix=None,
    )
    result = run_igvf_perturbseq(args)
    assert result["n_rows"] == 24

    signed = out_dir / "igvf_perturbseq_signed_term_gene.tsv"
    rows = list(csv.DictReader(signed.open("r", encoding="utf-8"), delimiter="\t"))
    assert {r["term"] for r in rows} == set(PERTS)
    # GENE_A_KO spikes TP53/MYC/EGFR up and KRAS/BRCA1/BRCA2 down.
    a_up = {r["gene_symbol"] for r in rows if r["term"] == "GENE_A_KO" and r["sign"] == "1"}
    a_dn = {r["gene_symbol"] for r in rows if r["term"] == "GENE_A_KO" and r["sign"] == "-1"}
    assert a_up == {"TP53", "MYC", "EGFR"}
    assert a_dn == {"KRAS", "BRCA1", "BRCA2"}

    graph_path = out_dir / "igvf_perturbseq_signed_term_gene.provenance_graph.json"
    assert graph_path.exists()
    payload = json.loads(graph_path.read_text(encoding="utf-8"))
    graph = next(iter(payload.values()))
    analyses = {n["id"].split(":")[1] for n in graph["nodes"] if n.get("type") == "AnalysisType"}
    assert "igvf_perturbseq" in analyses


def test_igvf_perturbseq_cli_entrypoint(tmp_path: Path):
    expr = tmp_path / "igvf_signatures.tsv"
    _write_perturbation_by_gene_matrix(expr)
    out_dir = tmp_path / "workflow_cli"
    code = main(
        [
            "workflows", "igvf_perturbseq",
            "--expression_tsv", str(expr),
            "--out_dir", str(out_dir),
            "--orientation", "perturbation_by_gene",
            "--z_threshold", "1.0",
            "--min_gmt_size", "2",
            "--gmt_name", "gene_set_library_crisp.gmt",
        ]
    )
    assert code == 0
    assert (out_dir / "igvf_perturbseq_signed_term_gene.tsv").exists()
    assert (out_dir / "gene_set_library_crisp.gmt").exists()


def test_igvf_perturbseq_long_de_preserves_legacy_filtering_and_ranking(tmp_path: Path):
    expression = tmp_path / "long_de.tsv"
    expression.write_text(
        "term\tgene\tgene_id\teffect\tpvalue\n"
        "P1\ta\tID1\t1\t0.01\n"
        "P1\ta\tID1\t3\t0.01\n"
        "P1\tb\tID2\t2\t0.01\n"
        "P1\tc\tID3\t-4\t0.01\n"
        "P1\td\tID4\t-1\t0.01\n"
        "P1\te\tID5\t9\t0.5\n",
        encoding="utf-8",
    )
    out_dir = tmp_path / "workflow_long_de"
    args = SimpleNamespace(
        expression_tsv=str(expression), mapping_file=None, out_dir=str(out_dir), organism="human", genome_build="hg38",
        orientation="perturbation_by_gene", gmt_name="gene_set_library_crisp.gmt", z_threshold=3.0, min_gmt_size=1,
        input_mode="long_de", sep="auto", term_column="term", gene_symbol_column="gene", gene_id_column="gene_id",
        effect_column="effect", ratio_column=None, score_column=None, pvalue_column="pvalue", pvalue_max=0.05,
        score_threshold=None, top_k_per_direction=1, provenance_overlay_json=None,
        provenance_mirror_local_prefix=None, provenance_mirror_remote_prefix=None,
    )
    assert run_igvf_perturbseq(args)["n_rows"] == 2
    signed = list(csv.DictReader((out_dir / "igvf_perturbseq_signed_term_gene.tsv").open(encoding="utf-8"), delimiter="\t"))
    assert {(row["gene_symbol"], row["sign"]) for row in signed} == {("A", "1"), ("C", "-1")}
