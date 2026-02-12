from omics2geneset.core.models import GeneWeights


def test_geneweights_sort_desc_tie_breaks_by_gene_id():
    rows = [
        {"gene_id": "G2", "weight": 1.0},
        {"gene_id": "G1", "weight": 1.0},
        {"gene_id": "G3", "weight": 2.0},
    ]
    gw = GeneWeights(rows).sort_desc()
    assert [r["gene_id"] for r in gw.rows] == ["G3", "G1", "G2"]
