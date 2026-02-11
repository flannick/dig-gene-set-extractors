from pathlib import Path

from omics2geneset.core.peak_to_gene import link_distance_decay, link_nearest_tss, link_promoter_overlap
from omics2geneset.io.bed import read_bed
from omics2geneset.io.gtf import read_genes_from_gtf


def test_link_methods_basic():
    genes = read_genes_from_gtf(Path("tests/data/toy.gtf"))
    peaks = read_bed(Path("tests/data/toy_peaks.bed"))

    prom = link_promoter_overlap(peaks, genes, promoter_upstream_bp=300, promoter_downstream_bp=100)
    assert any(l["gene_id"] == "GENE1" for l in prom)

    nearest = link_nearest_tss(peaks, genes, max_distance_bp=10000)
    assert len(nearest) == len(peaks)

    decay = link_distance_decay(peaks, genes, max_distance_bp=10000, decay_length_bp=1000, max_genes_per_peak=2)
    p0 = [d for d in decay if d["peak_index"] == 0]
    assert p0
    if len(p0) > 1:
        assert p0[0]["link_weight"] >= p0[-1]["link_weight"]
