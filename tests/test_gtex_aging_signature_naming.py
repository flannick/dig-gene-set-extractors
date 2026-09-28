from geneset_extractors.workflows.gtex_aging_signatures import _sanitize_tissue_label


def test_gtex_aging_signature_uses_underscore_separated_tissue_token() -> None:
    assert _sanitize_tissue_label("Adipose Tissue") == "Adipose_Tissue"
    assert _sanitize_tissue_label("Adipose - Subcutaneous") == "Adipose_Subcutaneous"
