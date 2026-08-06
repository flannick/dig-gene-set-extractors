# Testing DIG submission contracts

Run `pytest -q tests/test_submission_contract.py` after installing
`python -m pip install -e '.[dev]'`. The `rna_deg` contract uses the packaged
tiny `geneset_extractors/resources/submission_toy_deg.tsv` fixture and verifies the standard final files:
`geneset.tsv`, `geneset.meta.json`, and `geneset.provenance.json`.

Other registered entries are checked for CLI registration and module
importability. They intentionally do not all run in this low-cost interface:
some workflows require assay-specific input structures or optional runtimes.
