# ATAC Reference Bundle Setup

This guide covers how to use a versioned ATAC reference bundle with `omics2geneset`.
The primary workflow is direct bundle usage (no separate resource fetch step).

## 1) Download and unpack the bundle

If you have a published URL:

```bash
curl -fsSL -o /tmp/dig-gene-set-extractors-atac-refdata-v1.1.0.tar.gz <PUBLISHED_TARBALL_URL>
mkdir -p /tmp/dig-atac-refdata
tar -xzf /tmp/dig-gene-set-extractors-atac-refdata-v1.1.0.tar.gz -C /tmp/dig-atac-refdata
```

If you built it locally in this repo:

```bash
mkdir -p /tmp/dig-atac-refdata
tar -xzf refdata_build/dig-gene-set-extractors-atac-refdata-v1.1.0.tar.gz -C /tmp/dig-atac-refdata
```

After unpacking, your bundle root should be:

`/tmp/dig-atac-refdata/bundle`

Human builds supported by this bundle:

- `hg19` (`GRCh37`)
- `hg38` (`GRCh38`)

## 2) Create a local resources manifest (direct mode)

This helper writes a manifest where `filename` points inside the extracted bundle tree.
No additional downloads are required after this.

```bash
python scripts/make_local_resources_manifest.py \
  --bundle-root /tmp/dig-atac-refdata/bundle \
  --layout direct \
  --out /tmp/omics2geneset.local_resources.json
```

## 3) Run directly against the extracted bundle directory

```bash
omics2geneset resources status \
  --manifest /tmp/omics2geneset.local_resources.json \
  --manifest_mode replace \
  --resources_dir /tmp/dig-atac-refdata/bundle \
  --check_schema \
  --fast
```

```bash
omics2geneset convert atac_bulk \
  --peaks <peaks.bed> \
  --gtf <genes.gtf> \
  --out_dir <out_dir> \
  --organism human \
  --genome_build hg19 \
  --program_preset default \
  --resource_policy fail \
  --resources_manifest /tmp/omics2geneset.local_resources.json \
  --resources_dir /tmp/dig-atac-refdata/bundle
```

With `--use_reference_bundle true` (default), resource-backed methods are included automatically
unless `--program_methods` is explicitly set.
The converter auto-selects `*_hg19` or `*_hg38` resource IDs from `--genome_build`.

Optional: set the resource root once via environment variable:

```bash
export OMICS2GENESET_RESOURCES_DIR=/tmp/dig-atac-refdata/bundle
```

Then you can omit `--resources_dir` in convert/status commands.

## 4) Optional cache staging workflow

If you prefer to copy resources into the standard cache, generate a cache-style manifest and fetch:

```bash
python scripts/make_local_resources_manifest.py \
  --bundle-root /tmp/dig-atac-refdata/bundle \
  --layout cache \
  --out /tmp/omics2geneset.local_resources.cache.json

omics2geneset resources fetch \
  --manifest /tmp/omics2geneset.local_resources.cache.json \
  --manifest_mode replace \
  --preset atac_default_optional_hg19

omics2geneset resources fetch \
  --manifest /tmp/omics2geneset.local_resources.cache.json \
  --manifest_mode replace \
  --preset atac_default_optional_hg38
```

## Notes

- With direct mode, the only URL download is the initial bundle `.tar.gz`.
- `resources fetch` accepts HTTP(S), `file://`, and plain Unix paths in manifest `url` fields.
- `--manifest_mode replace` uses only your local manifest; `overlay` merges with bundled entries.
- Recommended for production runs: `--resource_policy fail` after `resources status --check_schema` is clean.
