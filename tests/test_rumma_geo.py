from __future__ import annotations

import json
from argparse import Namespace
from pathlib import Path

from geneset_extractors.extractors.converters.rumma_geo import run


def _write(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def test_reconstructs_notebook_filtering_and_duplicate_policy(tmp_path: Path) -> None:
    human_gmt, mouse_gmt = tmp_path / "human.gmt", tmp_path / "mouse.gmt"
    _write(human_gmt, "GSE1-0-vs-1-human up\t\tA\tA\tB\tC\tD\tE\nGSE1-0-vs-1-human dn\t\tF\tG\tH\tI\tJ\n")
    _write(mouse_gmt, "GSE2-0-vs-1-mouse dn\t\tM1\tM2\tM3\tM4\tM5\n")
    _write(tmp_path / "human.tsv", "#tax_id\tGeneID\tSymbol\ttype_of_gene\n9606\t1\tA\tprotein-coding\n9606\t2\tB\tprotein-coding\n9606\t3\tC\tprotein-coding\n9606\t4\tD\tprotein-coding\n9606\t5\tE\tprotein-coding\n9606\t6\tF\tprotein-coding\n9606\t7\tG\tprotein-coding\n9606\t8\tH\tprotein-coding\n9606\t9\tI\tprotein-coding\n9606\t10\tJ\tprotein-coding\n9606\t11\tK\tprotein-coding\n9606\t12\tL\tprotein-coding\n")
    _write(tmp_path / "mouse.tsv", "#tax_id\tGeneID\tSymbol\n10090\t101\tM1\n10090\t102\tM2\n10090\t103\tM3\n10090\t104\tM4\n10090\t105\tM5\n")
    _write(tmp_path / "orth.tsv", "#tax_id\tGeneID\tOther_tax_id\tOther_GeneID\n9606\t8\t10090\t101\n9606\t9\t10090\t102\n9606\t10\t10090\t103\n9606\t11\t10090\t104\n9606\t12\t10090\t105\n")
    _write(tmp_path / "selection.tsv", "source_term\tmodel_id\tstatus\tnormalized_term\tsearch_term\nGSE1-0-vs-1-human\tHZ2\tsignature\tGSE1_KO_human\tKO\nGSE2-0-vs-1-mouse\tHZ2\treversed\tGSE2_KO_mouse\tKO\n")
    roles = ["human_rummageo_gmt", "mouse_rummageo_gmt", "recorded_selection_manifest", "ncbi_human_gene_info", "ncbi_mouse_gene_info", "ncbi_gene_orthologs"]
    (tmp_path / "sources.json").write_text(json.dumps({"sources": {role: {"url": f"https://example.org/{role}", "version": "fixture-v1"} for role in roles}}), encoding="utf-8")
    result = run(Namespace(human_gmt=str(human_gmt), mouse_gmt=str(mouse_gmt), selection_manifest=str(tmp_path / "selection.tsv"), source_manifest=str(tmp_path / "sources.json"), human_gene_info=str(tmp_path / "human.tsv"), mouse_gene_info=str(tmp_path / "mouse.tsv"), gene_orthologs=str(tmp_path / "orth.tsv"), model_id="HZ2", out_dir=str(tmp_path / "out"), genome_build="hg38", min_genes=5, gmt_description="test", legacy_gmt=None, provenance_overlay_json=None, provenance_mirror_local_prefix=None, provenance_mirror_remote_prefix=None))
    assert result["n_sets"] == 2
    lines = (tmp_path / "out/genesets.gmt").read_text(encoding="utf-8").splitlines()
    assert [line.split("\t")[0] for line in lines] == ["GSE1_KO_human_dn", "GSE2_KO_mouse_up"]
    assert "A" not in lines[0]  # duplicate term/gene memberships are removed, not collapsed.
    assert lines[1].split("\t")[2:] == ["H", "I", "J", "K", "L"]
    metadata = json.loads((tmp_path / "out/geneset.meta.json").read_text(encoding="utf-8"))
    assert metadata["input"]["files"][0]["version"] == "fixture-v1"
