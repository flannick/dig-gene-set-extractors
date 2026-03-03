from geneset_extractors.cli import build_parser as build_parser_alias
from geneset_extractors.cli import build_parser as build_parser_legacy


def test_neutral_namespace_parser_alias() -> None:
    assert build_parser_alias is build_parser_legacy


def test_neutral_namespace_module_alias() -> None:
    from geneset_extractors.core.validate import validate_output_dir as alias_validate
    from geneset_extractors.core.validate import validate_output_dir as legacy_validate

    assert alias_validate is legacy_validate

