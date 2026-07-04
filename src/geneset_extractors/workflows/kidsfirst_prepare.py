#!/usr/bin/env python3
"""KidsFirst input-preparation workflow (DIG-owned).

Formalizes KidsFirst / CBTN tumor + normal RNA-seq matrix preparation into a single
reusable DIG workflow: build the tumor RSEM count matrix, extract the GTEx normal
tissue matrix (or accept pre-built matrices), align gene IDs across the two, and emit
the ``combined_counts.tsv`` + ``sample_metadata.tsv`` that feed
``geneset-extractors workflows rna_de_prepare``.

This logic was migrated verbatim from
geneset-extractor-dev/KidsFirst/src/{build_rsem_matrix,extract_gtex_counts,prepare_de_inputs}.py
so that DIG owns the reusable KidsFirst workflow and geneset-extractor-dev/KidsFirst
acts as the config/wrapper layer (branch two-repo standard).

Invoke via the DIG CLI:
    geneset-extractors workflows kidsfirst_prepare \
        --rsem_dir <dir> --manifest_tsv <tsv> \
        --gtex_gct <gct.gz> --gtex_sample_attrs <attrs.txt> --gtex_tissue "Whole Blood" \
        --study_id KF-TALL --out_dir <out>
or supply pre-built matrices with --tumor_counts / --normal_counts.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


# ---------------------------------------------------------------------------
# Tumor RSEM count matrix  (from build_rsem_matrix.py)
# ---------------------------------------------------------------------------
def _load_manifest(
    manifest_path: Path,
    filter_column: str | None = None,
    filter_value: str | None = None,
) -> dict[str, str]:
    """Returns {file_name: sample_id}. Optionally filter rows (e.g. CBTN diagnosis)."""
    mapping: dict[str, str] = {}
    skipped = 0
    with open(manifest_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if filter_column and filter_value:
                if row.get(filter_column, "").strip() != filter_value:
                    skipped += 1
                    continue
            fname = (row.get("file_name") or "").strip()
            sid = (row.get("sample_id") or "").strip()
            if fname and sid:
                mapping[fname] = sid
    if filter_column and filter_value:
        print(f"  Manifest filter: {filter_column}=={filter_value!r} "
              f"→ {len(mapping)} kept, {skipped} skipped", file=sys.stderr)
    return mapping


def _read_one_rsem(task: tuple[Path, str]) -> tuple[str, list[str], list[float]]:
    path, sample_id = task
    gene_ids: list[str] = []
    counts: list[float] = []
    open_fn = gzip.open if path.suffix == ".gz" else open
    with open_fn(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene_ids.append(row["gene_id"])
            counts.append(float(row["expected_count"]))
    return sample_id, gene_ids, counts


def build_tumor_matrix(
    rsem_dir: Path,
    manifest_path: Path | None,
    out_path: Path,
    workers: int = 4,
    filter_column: str | None = None,
    filter_value: str | None = None,
) -> None:
    rsem_files = sorted(rsem_dir.glob("*.rsem.genes.results.gz"))
    if not rsem_files:
        raise SystemExit(f"ERROR: no .rsem.genes.results.gz files in {rsem_dir}")

    file_to_sample: dict[str, str] = {}
    if manifest_path and manifest_path.exists():
        file_to_sample = _load_manifest(manifest_path, filter_column, filter_value)

    tasks: list[tuple[Path, str]] = []
    for f in rsem_files:
        if file_to_sample and f.name not in file_to_sample:
            continue
        sid = file_to_sample.get(f.name, f.name.split(".rsem.genes.results")[0])
        tasks.append((f, sid))

    print(f"Reading {len(tasks)} RSEM files with {workers} workers...", file=sys.stderr)
    results: dict[str, list[float]] = {}
    gene_ids: list[str] | None = None
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_read_one_rsem, t): t[1] for t in tasks}
        done = 0
        for future in as_completed(futures):
            sid, gids, counts = future.result()
            if gene_ids is None:
                gene_ids = gids
            results[sid] = counts
            done += 1
            if done % 100 == 0 or done == len(tasks):
                print(f"  {done}/{len(tasks)}", file=sys.stderr, end="\r")
    print(f"\n  {len(results)} samples, {len(gene_ids or [])} genes", file=sys.stderr)

    sample_ids = sorted(results.keys())
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["gene_id"] + sample_ids)
        for i, gene_id in enumerate(gene_ids or []):
            writer.writerow([gene_id] + [str(int(round(results[sid][i]))) for sid in sample_ids])
    print(f"Written tumor matrix: {out_path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# GTEx normal tissue matrix  (from extract_gtex_counts.py)
# ---------------------------------------------------------------------------
def _open(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def _get_tissue_sample_ids(sample_attrs_path: Path, tissue: str) -> set[str]:
    ids: set[str] = set()
    with open(sample_attrs_path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("SMTSD", "").strip() == tissue:
                sid = row.get("SAMPID", "").strip()
                if sid:
                    ids.add(sid)
    return ids


def extract_gtex_matrix(gct_path: Path, sample_attrs_path: Path, tissue: str, out_path: Path) -> None:
    tissue_ids = _get_tissue_sample_ids(sample_attrs_path, tissue)
    if not tissue_ids:
        raise SystemExit(f"ERROR: no samples found for tissue '{tissue}'")
    print(f"  {len(tissue_ids)} samples in SampleAttributes for '{tissue}'", file=sys.stderr)

    with _open(gct_path) as fh:
        fh.readline(); fh.readline()
        header = next(csv.reader(fh, delimiter="\t"))
    all_sample_ids = header[2:]
    keep_idx = [i for i, sid in enumerate(all_sample_ids) if sid in tissue_ids]
    keep_ids = [all_sample_ids[i] for i in keep_idx]
    if not keep_ids:
        raise SystemExit(f"ERROR: none of the {len(tissue_ids)} tissue samples found in GCT header")
    print(f"  {len(keep_ids)} samples found in GCT", file=sys.stderr)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_genes = 0
    with _open(gct_path) as fh, open(out_path, "w", newline="") as out_fh:
        fh.readline(); fh.readline()
        reader = csv.reader(fh, delimiter="\t")
        next(reader)
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow(["gene_id"] + keep_ids)
        for row in reader:
            if not row:
                continue
            counts = [row[2 + i] for i in keep_idx]
            writer.writerow([row[0].strip()] + counts)
            n_genes += 1
    print(f"  {n_genes} genes written\nWritten normal matrix: {out_path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Combined DE inputs  (from prepare_de_inputs.py)
# ---------------------------------------------------------------------------
def _strip_version(gene_id: str) -> str:
    return gene_id.split(".")[0] if "." in gene_id else gene_id


def _read_matrix(path: Path) -> tuple[list[str], list[str], dict[str, list[str]]]:
    gene_ids: list[str] = []
    data: dict[str, list[str]] = {}
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        sample_ids = header[1:]
        for row in reader:
            if not row:
                continue
            gid = _strip_version(row[0].strip())
            gene_ids.append(gid)
            data[gid] = row[1:]
    return gene_ids, sample_ids, data


def _read_tumor_metadata(path: Path, study_id: str) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sid = (row.get("Sample ID") or row.get("sample_id") or "").strip()
            diag = (row.get("Diagnosis (Source Text)") or row.get("diagnosis") or "").strip()
            if sid:
                mapping[sid] = diag or study_id
    return mapping


def _load_gene_map(path: Path | None) -> dict[str, str]:
    if path is None or not path.exists():
        return {}
    mapping: dict[str, str] = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gid = row.get("gene_id", "").strip()
            sym = row.get("gene_symbol", "").strip()
            if gid and sym:
                mapping[gid] = sym
    return mapping


def merge_de_inputs(
    tumor_counts_path: Path,
    normal_counts_path: Path,
    tumor_metadata_path: Path | None,
    study_id: str,
    out_dir: Path,
    normal_source: str = "GTEx",
    gene_map_path: Path | None = None,
) -> dict[str, int]:
    gene_map = _load_gene_map(gene_map_path)
    if gene_map:
        print(f"  Gene symbol map: {len(gene_map)} entries loaded", file=sys.stderr)

    tumor_genes, tumor_samples, tumor_data = _read_matrix(tumor_counts_path)
    print(f"  {len(tumor_samples)} tumor samples, {len(tumor_genes)} genes", file=sys.stderr)
    normal_genes, normal_samples, normal_data = _read_matrix(normal_counts_path)
    print(f"  {len(normal_samples)} normal samples, {len(normal_genes)} genes", file=sys.stderr)

    normal_gene_set = set(normal_genes)
    shared_genes = [g for g in tumor_genes if g in normal_gene_set]
    print(f"  {len(shared_genes)} shared genes after intersection", file=sys.stderr)
    if len(shared_genes) < 10000:
        print("WARNING: fewer than 10k shared genes — check gene ID format", file=sys.stderr)

    tumor_diag: dict[str, str] = {}
    if tumor_metadata_path and tumor_metadata_path.exists():
        tumor_diag = _read_tumor_metadata(tumor_metadata_path, study_id)

    all_samples = list(tumor_samples) + list(normal_samples)
    out_dir.mkdir(parents=True, exist_ok=True)

    counts_path = out_dir / "combined_counts.tsv"
    has_symbols = bool(gene_map)
    header = ["gene_id"] + (["gene_symbol"] if has_symbols else []) + all_samples
    with open(counts_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        for gene_id in shared_genes:
            tumor_row = tumor_data.get(gene_id, ["0"] * len(tumor_samples))
            normal_row = normal_data.get(gene_id, ["0"] * len(normal_samples))
            sym_cols = [gene_map.get(gene_id, "")] if has_symbols else []
            writer.writerow([gene_id] + sym_cols + tumor_row + normal_row)

    meta_path = out_dir / "sample_metadata.tsv"
    with open(meta_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t",
                                fieldnames=["sample_id", "condition", "source", "diagnosis"],
                                lineterminator="\n")
        writer.writeheader()
        for sid in tumor_samples:
            writer.writerow({"sample_id": sid, "condition": "tumor", "source": study_id,
                             "diagnosis": tumor_diag.get(sid, study_id)})
        for sid in normal_samples:
            writer.writerow({"sample_id": sid, "condition": "normal", "source": normal_source,
                             "diagnosis": "normal"})

    print(f"Done. combined_counts.tsv: {len(all_samples)} samples x {len(shared_genes)} genes; "
          f"sample_metadata.tsv: {len(tumor_samples)} tumor + {len(normal_samples)} normal", file=sys.stderr)
    return {"n_tumor": len(tumor_samples), "n_normal": len(normal_samples), "n_genes": len(shared_genes)}


# ---------------------------------------------------------------------------
# Workflow entry point
# ---------------------------------------------------------------------------
def run(args: argparse.Namespace) -> dict:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if getattr(args, "tumor_counts", None):
        tumor_counts = Path(args.tumor_counts)
    else:
        if not getattr(args, "rsem_dir", None):
            raise SystemExit("kidsfirst_prepare: provide --tumor_counts or --rsem_dir")
        tumor_counts = out_dir / "tumor_counts.tsv"
        build_tumor_matrix(
            Path(args.rsem_dir),
            Path(args.manifest_tsv) if args.manifest_tsv else None,
            tumor_counts, args.workers, args.filter_column, args.filter_value,
        )

    if getattr(args, "normal_counts", None):
        normal_counts = Path(args.normal_counts)
    else:
        if not (getattr(args, "gtex_gct", None) and getattr(args, "gtex_sample_attrs", None) and getattr(args, "gtex_tissue", None)):
            raise SystemExit("kidsfirst_prepare: provide --normal_counts or --gtex_gct/--gtex_sample_attrs/--gtex_tissue")
        normal_counts = out_dir / "normal_counts.tsv"
        extract_gtex_matrix(Path(args.gtex_gct), Path(args.gtex_sample_attrs), args.gtex_tissue, normal_counts)

    stats = merge_de_inputs(
        tumor_counts, normal_counts,
        Path(args.tumor_metadata) if args.tumor_metadata else None,
        args.study_id, out_dir,
        normal_source=(args.normal_source or "GTEx"),
        gene_map_path=Path(args.gene_map_tsv) if args.gene_map_tsv else None,
    )
    return {
        "out_dir": str(out_dir), "study_id": args.study_id,
        "combined_counts": str(out_dir / "combined_counts.tsv"),
        "sample_metadata": str(out_dir / "sample_metadata.tsv"),
        **stats,
    }


def add_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--study_id", required=True)
    # tumor source (build from RSEM, or supply a pre-built matrix)
    parser.add_argument("--rsem_dir", help="Directory of *.rsem.genes.results.gz")
    parser.add_argument("--manifest_tsv", default=None, help="RSEM manifest (file_name, sample_id[, filter col])")
    parser.add_argument("--filter_column", default=None, help="Manifest column to filter (e.g. CBTN diagnosis)")
    parser.add_argument("--filter_value", default=None)
    parser.add_argument("--tumor_counts", default=None, help="Pre-built tumor count matrix (skips RSEM build)")
    parser.add_argument("--workers", type=int, default=4)
    # normal source (extract from GTEx, or supply a pre-built matrix)
    parser.add_argument("--gtex_gct", default=None, help="GTEx GCT .gz")
    parser.add_argument("--gtex_sample_attrs", default=None, help="GTEx SampleAttributesDS.txt")
    parser.add_argument("--gtex_tissue", default=None, help="SMTSD value, e.g. 'Whole Blood'")
    parser.add_argument("--normal_counts", default=None, help="Pre-built normal count matrix (skips GTEx extract)")
    parser.add_argument("--normal_source", default="GTEx", help="Label for the normal cohort in sample_metadata")
    # shared
    parser.add_argument("--tumor_metadata", default=None, help="Tumor sample metadata (Sample ID, Diagnosis)")
    parser.add_argument("--gene_map_tsv", default=None, help="TSV with gene_id, gene_symbol columns")


def main() -> int:
    parser = argparse.ArgumentParser(description="KidsFirst tumor+normal DE-input preparation.")
    add_flags(parser)
    result = run(parser.parse_args())
    print(f"kidsfirst_prepare_completed out={result['out_dir']} "
          f"n_tumor={result['n_tumor']} n_normal={result['n_normal']} n_genes={result['n_genes']}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
