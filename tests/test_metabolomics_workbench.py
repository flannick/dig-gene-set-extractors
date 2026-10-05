from __future__ import annotations

from argparse import Namespace
from pathlib import Path

from geneset_extractors.extractors.converters.metabolomics_workbench import run


def test_hz1_deduplicates_filters_and_sorts(tmp_path: Path) -> None:
    edges = tmp_path / "edges.tsv"
    edges.write_text("Gene\tGene ID\tMetabolite\tMetabolite ID\tThreshold\nG5\t5\tMetabolite A\tC00001\t1\nG2\t2\tMetabolite A\tC00001\t1\nG1\t1\tMetabolite A\tC00001\t1\nG4\t4\tMetabolite A\tC00001\t1\nG3\t3\tMetabolite A\tC00001\t1\nG1\t1\tMetabolite A\tC00001\t1\nG1\t1\tMetabolite B\tC00002\t1\nG2\t2\tMetabolite B\tC00002\t1\nG3\t3\tMetabolite B\tC00002\t1\nG4\t4\tMetabolite B\tC00002\t1\n", encoding="utf-8")
    out = tmp_path / "out"
    result = run(Namespace(edges=str(edges), out_dir=str(out), model_id="HZ1", min_genes=5, genome_build="hg38", gmt_description="MW smoke", provenance_overlay_json=None))
    assert result["n_gene_sets"] == 1 and result["n_memberships"] == 5
    assert (out / "genesets.gmt").read_text(encoding="utf-8") == "Metabolite A\tMW smoke\tG1\tG2\tG3\tG4\tG5\n"
