# ATAC Reference Bundle Setup

This guide covers how to use versioned ATAC reference bundles with `omics2geneset`.
The primary workflow is direct bundle usage (no separate resource fetch step).
Recommended default: use a build-specific bundle (`hg19` or `hg38`).

## 1) Download and unpack a build-specific bundle

If you have published URLs:

```bash
# hg19
curl -fsSL -o /tmp/dig-gene-set-extractors-atac-refdata-hg19-v1.1.0.tar.gz <PUBLISHED_HG19_TARBALL_URL>
mkdir -p /tmp/dig-atac-refdata-hg19
tar -xzf /tmp/dig-gene-set-extractors-atac-refdata-hg19-v1.1.0.tar.gz -C /tmp/dig-atac-refdata-hg19

# hg38
curl -fsSL -o /tmp/dig-gene-set-extractors-atac-refdata-hg38-v1.1.0.tar.gz <PUBLISHED_HG38_TARBALL_URL>
mkdir -p /tmp/dig-atac-refdata-hg38
tar -xzf /tmp/dig-gene-set-extractors-atac-refdata-hg38-v1.1.0.tar.gz -C /tmp/dig-atac-refdata-hg38
```

If you built split bundles locally in this repo:

```bash
mkdir -p /tmp/dig-atac-refdata-hg19 /tmp/dig-atac-refdata-hg38
tar -xzf refdata_build/dig-gene-set-extractors-atac-refdata-hg19-v1.1.0.tar.gz -C /tmp/dig-atac-refdata-hg19
tar -xzf refdata_build/dig-gene-set-extractors-atac-refdata-hg38-v1.1.0.tar.gz -C /tmp/dig-atac-refdata-hg38
```

Bundle roots:

- `/tmp/dig-atac-refdata-hg19/bundle`
- `/tmp/dig-atac-refdata-hg38/bundle`

## 2) Create a local resources manifest (direct mode)

This helper writes a manifest where `filename` points inside the extracted bundle tree.
No additional downloads are required after this.

```bash
# hg19
python scripts/make_local_resources_manifest.py \
  --bundle-root /tmp/dig-atac-refdata-hg19/bundle \
  --layout direct \
  --genome-build hg19 \
  --out /tmp/omics2geneset.local_resources.hg19.json

# hg38
python scripts/make_local_resources_manifest.py \
  --bundle-root /tmp/dig-atac-refdata-hg38/bundle \
  --layout direct \
  --genome-build hg38 \
  --out /tmp/omics2geneset.local_resources.hg38.json
```

Use the manifest matching your `--genome_build` in `omics2geneset convert`.

## 3) Run directly against the extracted bundle directory

```bash
# hg19 example
omics2geneset resources status \
  --manifest /tmp/omics2geneset.local_resources.hg19.json \
  --manifest_mode replace \
  --resources_dir /tmp/dig-atac-refdata-hg19/bundle \
  --check_schema \
  --fast
```

```bash
# hg19 example
omics2geneset convert atac_bulk \
  --peaks <peaks.bed> \
  --gtf <genes.gtf> \
  --out_dir <out_dir> \
  --organism human \
  --genome_build hg19 \
  --program_preset default \
  --resource_policy fail \
  --resources_manifest /tmp/omics2geneset.local_resources.hg19.json \
  --resources_dir /tmp/dig-atac-refdata-hg19/bundle
```

With `--use_reference_bundle true` (default), resource-backed methods are included automatically unless `--program_methods` is explicitly set.
The converter auto-selects `*_hg19` or `*_hg38` resource IDs from `--genome_build`.

Optional environment variable:

```bash
export OMICS2GENESET_RESOURCES_DIR=/tmp/dig-atac-refdata-hg19/bundle
```

## 4) Optional cache staging workflow

If you prefer to copy resources into the standard cache, generate a cache-style manifest and fetch:

```bash
# hg19
python scripts/make_local_resources_manifest.py \
  --bundle-root /tmp/dig-atac-refdata-hg19/bundle \
  --layout cache \
  --genome-build hg19 \
  --out /tmp/omics2geneset.local_resources.hg19.cache.json

omics2geneset resources fetch \
  --manifest /tmp/omics2geneset.local_resources.hg19.cache.json \
  --manifest_mode replace \
  --preset atac_default_optional_hg19
```

```bash
# hg38
python scripts/make_local_resources_manifest.py \
  --bundle-root /tmp/dig-atac-refdata-hg38/bundle \
  --layout cache \
  --genome-build hg38 \
  --out /tmp/omics2geneset.local_resources.hg38.cache.json

omics2geneset resources fetch \
  --manifest /tmp/omics2geneset.local_resources.hg38.cache.json \
  --manifest_mode replace \
  --preset atac_default_optional_hg38
```

## 5) Create split bundles locally

If you already have a combined human bundle extracted at `refdata_build/bundle`:

```bash
python scripts/split_human_atac_reference_bundles.py \
  --source-bundle refdata_build/bundle \
  --out-dir refdata_build
```

This writes:

- `refdata_build/dig-gene-set-extractors-atac-refdata-hg19-v1.1.0.tar.gz`
- `refdata_build/dig-gene-set-extractors-atac-refdata-hg38-v1.1.0.tar.gz`
- `refdata_build/SHA256SUMS.human_split.txt`

## Notes

- With direct mode, the only URL download is the initial bundle `.tar.gz`.
- `resources fetch` accepts HTTP(S), `file://`, and plain Unix paths in manifest `url` fields.
- `--manifest_mode replace` uses only your local manifest; `overlay` merges with bundled entries.
- Recommended for production runs: `--resource_policy fail` after `resources status --check_schema` is clean.
