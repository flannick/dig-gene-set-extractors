from geneset_extractors.workflows.gtex_aging_signatures import _build_comparisons, _sanitize_tissue_label


def test_gtex_aging_signature_uses_underscore_separated_tissue_token() -> None:
    assert _sanitize_tissue_label("Adipose Tissue") == "Adipose_Tissue"
    assert _sanitize_tissue_label("Adipose - Subcutaneous") == "Adipose_Subcutaneous"


def test_gtex_aging_signature_emits_readable_age_pair_for_comparison_directory() -> None:
    comparisons, _selected_samples, _audit = _build_comparisons(
        [
            *({"sample_id": f"reference_{index}", "age_group": "20-29"} for index in range(3)),
            *({"sample_id": f"case_{index}", "age_group": "30-39"} for index in range(3)),
        ],
        tissue_label="Adipose Tissue",
        reference_age_group="20-29",
        comparison_age_groups=["30-39"],
        random_state=1,
        min_samples_per_group=3,
    )
    assert comparisons[0]["gmt_comparison_label"] == "20-29_30-39"
