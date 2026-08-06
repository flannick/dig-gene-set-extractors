# Workflow and converter contract

Substantive source-data preprocessing, statistical analysis, normalization,
differential testing, gene mapping, ranking, gene-set construction, and
reusable converters belong in `dig-gene-set-extractors`.

The `submission` CLI derives contracts from the existing converter registry and
the existing workflow CLI parser. Each JSON record reports identifier, source
module, CLI command, assay type, available input/output contract information,
fixture/smoke-test availability, and API stability. The wrapper repository may
configure and dispatch these interfaces but must not duplicate their logic.
