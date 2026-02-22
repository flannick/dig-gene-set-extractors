# ATAC Reference Bundle Setup

This guide covers how to use a versioned ATAC reference bundle with `omics2geneset`.

## 1) Download and unpack the bundle

If you have a published URL:

```bash
curl -fsSL -o /tmp/dig-gene-set-extractors-atac-refdata-v1.0.0.tar.gz <PUBLISHED_TARBALL_URL>
mkdir -p /tmp/dig-atac-refdata
tar -xzf /tmp/dig-gene-set-extractors-atac-refdata-v1.0.0.tar.gz -C /tmp/dig-atac-refdata
```

If you built it locally in this repo:

```bash
mkdir -p /tmp/dig-atac-refdata
tar -xzf refdata_build/dig-gene-set-extractors-atac-refdata-v1.0.0.tar.gz -C /tmp/dig-atac-refdata
```

After unpacking, your bundle root should be:

`/tmp/dig-atac-refdata/bundle`

## 2) Create a local resources manifest using Unix paths

This helper script writes a manifest where each resource `url` is a plain Unix file path (no URI required).

```bash
python scripts/make_local_resources_manifest.py \
  --bundle-root /tmp/dig-atac-refdata/bundle \
  --out /tmp/omics2geneset.local_resources.json
```

## 3) Fetch resources from those local paths into cache

```bash
omics2geneset resources fetch \
  --manifest /tmp/omics2geneset.local_resources.json \
  --manifest_mode replace \
  --preset atac_default_optional
```

Optional: specify an explicit cache directory.

```bash
omics2geneset resources fetch \
  --manifest /tmp/omics2geneset.local_resources.json \
  --manifest_mode replace \
  --resources_dir /tmp/omics2geneset-resources \
  --preset atac_default_optional
```

## 4) Use the manifest in conversion runs

```bash
omics2geneset convert atac_bulk \
  --peaks <peaks.bed> \
  --gtf <genes.gtf> \
  --out_dir <out_dir> \
  --organism human \
  --genome_build hg38 \
  --program_methods ref_ubiquity_penalty \
  --resources_manifest /tmp/omics2geneset.local_resources.json
```

## Notes

- `resources fetch` accepts HTTP(S), `file://`, and plain Unix paths in the resource `url` field.
- `--manifest_mode replace` uses only your local manifest; `overlay` merges with bundled entries.
- Verify local file availability with:

```bash
omics2geneset resources status \
  --manifest /tmp/omics2geneset.local_resources.json \
  --manifest_mode replace \
  --fast
```
