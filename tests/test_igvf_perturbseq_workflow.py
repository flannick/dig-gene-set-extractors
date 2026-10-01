from pathlib import Path
from types import SimpleNamespace

from geneset_extractors.workflows.igvf_perturbseq import run


def test_long_de_writes_signed_gmt(tmp_path: Path) -> None:
    source = tmp_path / "input.tsv"
    source.write_text("term\tgene\teffect\tp\nA\tG1\t2\t0.01\nA\tG2\t2\t0.01\nA\tG3\t2\t0.01\nA\tG4\t2\t0.01\nA\tG5\t2\t0.01\n", encoding="utf-8")
    args = SimpleNamespace(expression_tsv=str(source), out_dir=str(tmp_path / "out"), input_mode="long_de", sep="auto", term_column="term", gene_symbol_column="gene", gene_id_column=None, effect_column="effect", ratio_column=None, score_column=None, pvalue_column="p", pvalue_max=0.05, score_threshold=None, top_k_per_direction=200, gmt_name="out.gmt", min_gmt_size=5, z_threshold=3.0, orientation="perturbation_by_gene", mapping_file=None, provenance_overlay_json=None)
    assert run(args)["n_rows"] == 5
    assert (tmp_path / "out" / "out.gmt").is_file()
    assert (tmp_path / "out" / "igvf_perturbseq_signed_term_gene.provenance_graph.json").is_file()
