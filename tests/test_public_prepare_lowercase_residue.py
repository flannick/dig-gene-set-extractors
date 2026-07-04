from geneset_extractors.extractors.proteomics.public_prepare import (
    _parse_site_label,
    _assay_type_qc,
)


def test_parse_lowercase_single_site_normalizes_to_uppercase():
    parsed = _parse_site_label("NP_100001:s123", "phospho")
    assert parsed["site_parser_status"] == "single_site"
    assert parsed["residue"] == "S"
    assert parsed["position"] == "123"
    assert parsed["site_id"] == "NP_100001|S|123|phospho"


def test_parse_lowercase_multi_site_group():
    parsed = _parse_site_label("NP_100002:s10s14", "phospho")
    assert parsed["site_parser_status"] == "site_group"
    assert parsed["n_sites_in_group"] == 2


def test_assay_qc_counts_lowercase_sty_as_phospho_like():
    rows = [
        {"raw_site_label": "NP_1:s5"},
        {"raw_site_label": "NP_2:t9"},
        {"raw_site_label": "NP_3:y11"},
        {"raw_site_label": "NP_4:s20"},
    ]
    result = _assay_type_qc(
        rows,
        ptm_type="phospho",
        assay_type_policy="warn",
        min_phospho_like_fraction=0.6,
        max_k_fraction=0.25,
    )
    assert result["dominant_residue_family"] == "phospho_like"
    assert result["status"] == "pass"
