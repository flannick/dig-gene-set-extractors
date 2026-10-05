from __future__ import annotations

from argparse import Namespace
from pathlib import Path

from geneset_extractors.extractors.converters.impc import run


def test_impc_hz1_deduplicates_and_filters(tmp_path: Path) -> None:
    assertions = tmp_path / "assertions.csv"
    assertions.write_text("marker_symbol,mp_term_id,mp_term_name\nAdk,MP:0004882,enlarged lung\nCldn18,MP:0004882,enlarged lung\nKcnab3,MP:0004882,enlarged lung\nLrrc8a,MP:0004882,enlarged lung\nTwist2,MP:0004882,enlarged lung\nAdk,MP:0004882,enlarged lung\nPrrt2,MP:0001454,abnormal behavior\nKcne1,MP:0001454,abnormal behavior\nPtchd1,MP:0001454,abnormal behavior\nNab2,MP:0001454,abnormal behavior\n", encoding="utf-8")
    mapping = tmp_path / "mapping.tsv"
    mapping.write_text("ADK\tADK\nCLDN18\tCLDN18\nKCNAB3\tKCNAB3\nLRRC8A\tLRRC8A\nTWIST2\tTWIST2\nPRRT2\tPRRT2\nKCNE1\tKCNE1\nPTCHD1\tPTCHD1\nNAB2\tNAB2\n", encoding="utf-8")
    out = tmp_path / "out"
    result = run(Namespace(assertions=str(assertions), symbol_mapping=str(mapping), out_dir=str(out), model_id="HZ1", min_genes=5, genome_build="hg38", gmt_description="IMPC smoke", provenance_overlay_json=None))
    assert result["n_gene_sets"] == 1 and result["n_memberships"] == 5
    assert (out / "genesets.gmt").read_text(encoding="utf-8") == "Enlarged Lung (MP:0004882)\tIMPC smoke\tADK\tCLDN18\tKCNAB3\tLRRC8A\tTWIST2\n"
    assert (out / "geneset.meta.json").is_file()
