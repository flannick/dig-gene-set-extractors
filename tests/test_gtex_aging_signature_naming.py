from pathlib import Path


WORKFLOW_SOURCE = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "geneset_extractors"
    / "workflows"
    / "gtex_aging_signatures.py"
)


def test_gtex_aging_signature_source_uses_underscore_separated_tissue_token() -> None:
    source = WORKFLOW_SOURCE.read_text(encoding="utf-8")

    assert 'return "_".join(parts)' in source


def test_gtex_aging_signature_source_emits_comparison_label_for_gmt_layout() -> None:
    source = WORKFLOW_SOURCE.read_text(encoding="utf-8")

    assert 'gmt_comparison_label = f"{reference_age_group}_{age_group}"' in source
    assert '"gmt_comparison_label": gmt_comparison_label' in source
    assert '"comparison_id", "gmt_comparison_label", "comparison_kind"' in source
