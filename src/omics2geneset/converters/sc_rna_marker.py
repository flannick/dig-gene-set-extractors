from __future__ import annotations

import csv
from pathlib import Path

from omics2geneset.core.metadata import input_file_record, make_metadata, write_metadata
from omics2geneset.core.models import GeneWeights
from omics2geneset.core.normalization import normalize


def _read_counts(path: str | Path, gene_col: str, barcode_col: str, value_col: str) -> list[tuple[str, str, float]]:
    rows: list[tuple[str, str, float]] = []
    with Path(path).open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or gene_col not in reader.fieldnames or barcode_col not in reader.fieldnames or value_col not in reader.fieldnames:
            raise ValueError(f"counts_tsv must contain {gene_col}, {barcode_col}, {value_col}")
        for row in reader:
            rows.append((str(row[gene_col]), str(row[barcode_col]), float(row[value_col])))
    return rows


def _read_groups(path: str | Path, barcode_col: str, group_col: str) -> dict[str, str]:
    out: dict[str, str] = {}
    with Path(path).open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or barcode_col not in reader.fieldnames or group_col not in reader.fieldnames:
            raise ValueError(f"groups_tsv must contain {barcode_col}, {group_col}")
        for row in reader:
            out[str(row[barcode_col])] = str(row[group_col])
    return out


def _safe_group_name(name: str) -> str:
    return name.replace("/", "_").replace(" ", "_")


def _summarize(values: list[float], method: str, n_cells: int) -> float:
    if method == "sum_counts":
        return float(sum(values))
    if method == "mean_counts":
        return float(sum(values) / max(n_cells, 1))
    if method == "frac_cells_nonzero":
        return float(sum(1 for v in values if v > 0) / max(n_cells, 1))
    raise ValueError(f"Unsupported summary method: {method}")


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = _read_counts(args.counts_tsv, args.gene_id_column, args.barcode_column, args.value_column)
    barcodes = sorted({bc for _, bc, _ in rows})

    if args.groups_tsv:
        groups_map = _read_groups(args.groups_tsv, args.group_barcode_column, args.group_column)
        group_to_barcodes: dict[str, set[str]] = {}
        for bc in barcodes:
            group = groups_map.get(bc)
            if group is not None:
                group_to_barcodes.setdefault(group, set()).add(bc)
    else:
        group_to_barcodes = {"all": set(barcodes)}

    files = [input_file_record(args.counts_tsv, "counts_tsv")]
    if args.groups_tsv:
        files.append(input_file_record(args.groups_tsv, "groups_tsv"))

    manifest_rows: list[tuple[str, str]] = []
    for group, group_barcodes in group_to_barcodes.items():
        gene_to_values: dict[str, list[float]] = {}
        for gene, bc, val in rows:
            if bc in group_barcodes:
                gene_to_values.setdefault(gene, []).append(val)

        summarized = {gene: _summarize(vals, args.peak_summary, len(group_barcodes)) for gene, vals in gene_to_values.items()}
        final = normalize(summarized, args.normalize)
        geneset_rows = [{"gene_id": gid, "weight": w} for gid, w in final.items()]
        GeneWeights(geneset_rows).sort_desc().to_tsv(
            (out_dir / f"group={_safe_group_name(group)}" / "geneset.tsv") if args.groups_tsv else (out_dir / "geneset.tsv")
        )

        group_dir = out_dir / f"group={_safe_group_name(group)}" if args.groups_tsv else out_dir
        group_dir.mkdir(parents=True, exist_ok=True)

        meta = make_metadata(
            converter_name="sc_rna_marker",
            parameters={**vars(args).copy(), "group": group},
            data_type="rna_seq",
            assay="single_cell",
            organism=args.organism,
            genome_build=args.genome_build,
            files=files,
            gene_annotation={"gtf_path": "", "source": "none", "gene_id_field": args.gene_id_column},
            weights={
                "weight_type": "nonnegative",
                "normalization": {"method": args.normalize, "target_sum": 1.0 if args.normalize == "l1" else None},
                "aggregation": "sum",
            },
            summary={
                "n_input_features": len(rows),
                "n_genes": len(geneset_rows),
                "n_features_assigned": len(rows),
                "fraction_features_assigned": 1.0 if rows else 0.0,
            },
        )
        write_metadata(group_dir / "geneset.meta.json", meta)
        manifest_rows.append((group, str(group_dir)))

    if args.groups_tsv:
        with (out_dir / "manifest.tsv").open("w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["group", "path"])
            writer.writerows(manifest_rows)

    return {"n_peaks": len(rows), "n_genes": len({g for g, _, _ in rows}), "out_dir": str(out_dir), "n_groups": len(group_to_barcodes)}
