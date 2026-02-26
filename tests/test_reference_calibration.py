from omics2geneset.core.reference_calibration import peak_overlap_mask


def test_peak_overlap_mask_reports_overlaps_per_peak():
    peaks = [
        {"chrom": "chr1", "start": 100, "end": 200},
        {"chrom": "chr1", "start": 300, "end": 350},
        {"chrom": "chr2", "start": 50, "end": 80},
    ]
    ref = [
        {"chrom": "chr1", "start": 150, "end": 160, "idf_ref": 2.0},
        {"chrom": "chr2", "start": 1, "end": 10, "idf_ref": 3.0},
    ]

    assert peak_overlap_mask(peaks, ref) == [True, False, False]
