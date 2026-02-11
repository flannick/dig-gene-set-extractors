from __future__ import annotations

import csv
from pathlib import Path

from omics2geneset.core.metadata import input_file_record, make_metadata, write_metadata
from omics2geneset.core.models import GeneWeights
from omics2geneset.core.normalization import normalize
from omics2geneset.core.scoring import transform_peak_weight


def _read_scores(path: str | Path, gene_id_column: str, score_column: str) -> dict[str, float]:
    out: dict[str, float] = {}
    with Path(path).open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or gene_id_column not in reader.fieldnames or score_column not in reader.fieldnames:
            raise ValueError(f"Input TSV must contain {gene_id_column} and {score_column}")
        for row in reader:
            gid = str(row[gene_id_column])
            score = float(row[score_column])
            out[gid] = out.get(gid, 0.0) + score
    return out


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    raw = _read_scores(args.proteomics_tsv, args.gene_id_column, args.score_column)
    transformed = {gid: transform_peak_weight(val, args.score_transform) for gid, val in raw.items()}
    final = normalize(transformed, args.normalize)

    rows = [{"gene_id": gid, "weight": w} for gid, w in final.items()]
    GeneWeights(rows).sort_desc().to_tsv(out_dir / "geneset.tsv")

    meta = make_metadata(
        converter_name="proteomics_diff",
        parameters=vars(args).copy(),
        data_type="proteomics",
        assay="bulk",
        organism=args.organism,
        genome_build=args.genome_build,
        files=[input_file_record(args.proteomics_tsv, "proteomics_tsv")],
        gene_annotation={"gtf_path": "", "source": "none", "gene_id_field": args.gene_id_column},
        weights={
            "weight_type": "signed" if args.score_transform == "signed" else "nonnegative",
            "normalization": {"method": args.normalize, "target_sum": 1.0 if args.normalize == "l1" else None},
            "aggregation": "sum",
        },
        summary={
            "n_input_features": len(raw),
            "n_genes": len(rows),
            "n_features_assigned": len(raw),
            "fraction_features_assigned": 1.0 if raw else 0.0,
        },
    )
    write_metadata(out_dir / "geneset.meta.json", meta)
    return {"n_peaks": len(raw), "n_genes": len(rows), "out_dir": str(out_dir)}
