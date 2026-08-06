# DIG submission contract

This repository exposes a small, machine-readable view of its existing CLI
contracts for coordinated new-library reviews:

```bash
geneset-extractors submission list
geneset-extractors submission describe rna_deg
geneset-extractors submission validate rna_deg
```

It is not wrapper orchestration. See [workflow-contract.md](workflow-contract.md),
[testing.md](testing.md), and [paired-prs.md](paired-prs.md).

GitHub Actions runs this interface in the stable branch-protection check named
**`test-dig-submission-interface`**. Reproduce it locally with
`pytest -q tests/test_submission_contract.py`, then run the three commands
above from an environment with DIG installed.
