from pathlib import Path

from omics2geneset.core.gmt import choose_gene_tokens, sanitize_gmt_name, write_gmt


def test_choose_gene_tokens_deduplicates_with_symbol_preference():
    rows = [
        {"gene_id": "GENE1", "gene_symbol": "A", "score": 4.0},
        {"gene_id": "GENE2", "gene_symbol": "A", "score": 3.0},
        {"gene_id": "GENE3", "gene_symbol": "", "score": 2.0},
        {"gene_id": "GENE4", "score": 1.0},
    ]
    tokens = choose_gene_tokens(rows, prefer_symbol=True)
    assert tokens == ["A", "GENE3", "GENE4"]


def test_write_gmt_uses_single_tab_and_space_separated_genes(tmp_path: Path):
    out = tmp_path / "genesets.gmt"
    set_name = sanitize_gmt_name("my set / alpha")
    write_gmt([(set_name, ["G1", "G2", "G3"])], out)

    text = out.read_text(encoding="utf-8")
    lines = text.splitlines()
    assert len(lines) == 1
    assert lines[0].count("\t") == 1
    name, genes = lines[0].split("\t")
    assert name == "my_set___alpha"
    assert genes == "G1 G2 G3"
