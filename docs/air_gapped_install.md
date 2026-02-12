# Offline / Air-Gapped Install

This guide is for environments that cannot access PyPI.

## Common Online Case (for reference)

If you have internet access, use the normal install path instead:

```bash
python -m pip install -U pip
python -m pip install -e ".[dev]"
```

## Offline Setup Requirements

Prepare a local wheelhouse that includes at least these direct packages:

- `setuptools`
- `wheel`
- `build`
- `jsonschema`
- `pytest` (if installing with `.[dev]`)

Reference template:

- `offline/wheelhouse_manifest.example.txt`

## Recommended Offline Install Flow

```bash
scripts/check_wheelhouse.py --wheelhouse "$WHEELHOUSE" --with-dev
scripts/install_offline.sh --wheelhouse "$WHEELHOUSE" --python python
```

What this does:

1. Validates required direct wheels exist.
2. Bootstraps build prerequisites from the local wheelhouse.
3. Installs the package with `--no-build-isolation`.
4. Verifies CLI startup (`python -m omics2geneset.cli list`).

## Manual Offline Commands

```bash
python -m pip install --no-index --find-links "$WHEELHOUSE" setuptools wheel build
python -m pip install --no-index --find-links "$WHEELHOUSE" -e ".[dev]" --no-build-isolation
```

## Notes on `--no-deps`

`pip install -e . --no-build-isolation --no-deps` can fail in a fresh venv with:

- `invalid command 'bdist_wheel'`

Use this bootstrap first:

```bash
scripts/bootstrap_build_tools.sh --python python --wheelhouse "$WHEELHOUSE"
python -m pip install -e . --no-build-isolation --no-deps
```
