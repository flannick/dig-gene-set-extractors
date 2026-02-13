from pathlib import Path

from omics2geneset.core.gmt import build_gmt_sets_from_rows, choose_gene_tokens, sanitize_gmt_name, write_gmt


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


def test_choose_gene_tokens_require_symbol_drops_ensembl_like_symbols():
    rows = [
        {"gene_id": "ENSG000001.1", "gene_symbol": "ENSG000001", "score": 5.0},
        {"gene_id": "ENSG000002.1", "gene_symbol": "CXCR4", "score": 4.0},
    ]
    tokens = choose_gene_tokens(rows, prefer_symbol=True, require_symbol=True)
    assert tokens == ["CXCR4"]


def test_build_gmt_sets_filters_by_biotype_allowlist():
    rows = [
        {"gene_id": "ENSGA.1", "gene_symbol": "GENEA", "gene_biotype": "protein_coding", "score": 5.0},
        {"gene_id": "ENSGB.1", "gene_symbol": "GENEB", "gene_biotype": "lncRNA", "score": 4.0},
        {"gene_id": "ENSGC.1", "gene_symbol": "GENEC", "gene_biotype": "processed_pseudogene", "score": 3.0},
    ]
    sets, _ = build_gmt_sets_from_rows(
        rows=rows,
        base_name="demo",
        prefer_symbol=True,
        min_genes=1,
        max_genes=10,
        topk_list=[10],
        mass_list=[],
        split_signed=False,
        require_symbol=True,
        allowed_biotypes={"protein_coding", "lncrna"},
    )
    assert sets
    _, genes = sets[0]
    assert genes == ["GENEA", "GENEB"]
