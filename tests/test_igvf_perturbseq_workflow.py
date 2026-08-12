from __future__ import annotations

import csv
import gzip
import importlib
import sys
from pathlib import Path

from geneset_extractors.cli import main


def test_igvf_workflow_contract_import_does_not_require_pandas(monkeypatch) -> None:
    real_import = __import__

    def without_pandas(name, *args, **kwargs):
        if name == "pandas":
            raise ModuleNotFoundError("No module named 'pandas'")
        return real_import(name, *args, **kwargs)

    monkeypatch.delitem(sys.modules, "geneset_extractors.workflows.igvf_perturbseq", raising=False)
    monkeypatch.setattr("builtins.__import__", without_pandas)
    module = importlib.import_module("geneset_extractors.workflows.igvf_perturbseq")
    assert callable(module.run)


def test_igvf_perturbseq_preserves_signed_effect_and_top_k(tmp_path: Path) -> None:
    source = tmp_path / "input.tsv.gz"
    with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["term", "symbol", "gene_id", "effect", "p_value"], delimiter="\t")
        writer.writeheader()
        writer.writerows([
            {"term": "guide", "symbol": "A", "gene_id": "ENSG1", "effect": "2", "p_value": "0.01"},
            {"term": "guide", "symbol": "B", "gene_id": "ENSG2", "effect": "-3", "p_value": "0.01"},
            {"term": "guide", "symbol": "C", "gene_id": "ENSG3", "effect": "1", "p_value": "0.01"},
            {"term": "guide", "symbol": "DROP", "gene_id": "ENSG4", "effect": "9", "p_value": "0.2"},
        ])
    manifest = tmp_path / "analysis.tsv"
    manifest.write_text(
        "analysis_set_id\tsep\tterm_column\tgene_symbol_column\tgene_id_column\teffect_column\tratio_column\tscore_column\tpvalue_column\tpvalue_max\ttop_k_per_direction\n"
        "IGVFTEST\tauto\tterm\tsymbol\tgene_id\teffect\tNA\tNA\tp_value\t0.05\t1\n",
        encoding="utf-8",
    )
    out = tmp_path / "out"
    assert main(["workflows", "igvf_perturbseq", "--expression_tsv", str(source), "--analysis_set_manifest", str(manifest), "--analysis_set_id", "IGVFTEST", "--out_dir", str(out), "--min_gmt_size", "1"]) == 0
    signed = (out / "workflow" / "igvf_perturbseq_signed_term_gene.tsv").read_text(encoding="utf-8")
    assert "\tB\t3.0\t-1" in signed
    assert "DROP" not in signed
    gmt = (out / "extractor" / "genesets.gmt").read_text(encoding="utf-8")
    assert "IGVF_Perturb_Seq_guide_dn" in gmt
