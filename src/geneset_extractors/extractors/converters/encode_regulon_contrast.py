from __future__ import annotations

import collections
import json
import zipfile
from pathlib import Path

from geneset_extractors.core.gmt import write_gmt
from geneset_extractors.core.metadata import input_file_record, make_metadata, write_metadata
from geneset_extractors.core.provenance import activate_runtime_context

ENCODE_CITATION = (
    "ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the "
    "human genome. Nature. 2012;489(7414):57-74. doi:10.1038/nature11247. "
    "ENCODE/NHGRI; public data. https://www.encodeproject.org"
)

_Record = tuple[str, frozenset[str], str]
# (setdir_name, genes, target_label)


def _safe_name(s: str, max_len: int = 60) -> str:
    return "".join(c if c.isalnum() else "_" for c in s)[:max_len]


def _load_records_from_dir(src: Path) -> list[_Record]:
    recs: list[_Record] = []
    for sub in sorted(src.iterdir()):
        tsv = sub / "geneset.tsv"
        if not tsv.is_file():
            continue
        genes: set[str] = set()
        with tsv.open("r", encoding="utf-8") as fh:
            next(fh, None)
            for line in fh:
                g = line.strip()
                if g:
                    genes.add(g)
        if not genes:
            continue
        target = sub.name
        meta_path = sub / "geneset.meta.json"
        if meta_path.is_file():
            try:
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
                target = meta.get("target", target) or target
            except Exception:
                pass
        recs.append((sub.name, frozenset(genes), target))
    return recs


def _load_records_from_zip(zpath: Path) -> list[_Record]:
    recs: list[_Record] = []
    with zipfile.ZipFile(zpath, "r") as zf:
        names = zf.namelist()
        for n in names:
            if not n.endswith("geneset.tsv"):
                continue
            setdir = n.split("/")[-2]
            raw = zf.read(n).decode("utf-8", "replace").splitlines()
            genes: set[str] = {ln.strip() for ln in raw[1:] if ln.strip()}
            if not genes:
                continue
            target = setdir
            meta_name = n.rsplit("/", 1)[0] + "/geneset.meta.json"
            if meta_name in names:
                try:
                    meta = json.loads(zf.read(meta_name).decode("utf-8", "replace"))
                    target = meta.get("target", target) or target
                except Exception:
                    pass
            recs.append((setdir, frozenset(genes), target))
    return recs


def _write_contrast_set(
    set_dir: Path,
    set_name: str,
    description: str,
    genes: list[str],
    assay: str,
    target: str,
    direction: str,
    low_threshold: float,
    high_threshold: float,
    n_factors: int,
    input_file_rec: dict | None,
    regulon_src_dir: str | None = None,
) -> None:
    set_dir.mkdir(parents=True, exist_ok=True)
    with (set_dir / "geneset.tsv").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("gene\n")
        for g in genes:
            fh.write(f"{g}\n")
    write_gmt([(set_name, genes)], set_dir / "genesets.gmt")

    params: dict = {
        "assay": assay,
        "target": target,
        "direction": direction,
        "low_prevalence_threshold": low_threshold,
        "high_prevalence_threshold": high_threshold,
        "background_method": "leave_one_out_binding_prevalence",
        "background_n_factors": n_factors,
        "encode_citation": ENCODE_CITATION,
        "caveat": (
            "Cross-factor binding-prevalence background (relative binding specificity); "
            "removes promiscuous / HOT-region genes that are bound by many TFs/RBPs. "
            "Not a read-count differential (no DESeq2/edgeR)."
        ),
    }
    if regulon_src_dir is not None:
        params["regulon_src_dir"] = regulon_src_dir
    meta = make_metadata(
        converter_name="encode_regulon_contrast",
        parameters=params,
        data_type="transcription_factor_binding" if "ChIP" in assay else "rna_binding",
        assay=assay,
        organism="human",
        genome_build="GRCh38",
        files=[input_file_rec] if input_file_rec is not None else [],
        gene_annotation={
            "mode": "inherited",
            "source": "encode_regulon_output",
        },
        weights={
            "weight_type": "nonnegative",
            "normalization": {"method": "none", "target_sum": None},
            "aggregation": "loo_prevalence_threshold",
        },
        summary={
            "n_input_features": n_factors,
            "n_genes": len(genes),
            "n_features_assigned": len(genes),
            "fraction_features_assigned": 1.0,
            "n_sets_emitted": 1,
        },
        gene_set_description=description,
        output_files=[
            {"path": "genesets.gmt", "role": "gmt_library"},
            {"path": "geneset.tsv", "role": "selected_program"},
            {"path": "geneset.meta.json", "role": "metadata_json"},
        ],
    )
    write_metadata(set_dir / "geneset.meta.json", meta)


def run(args) -> dict[str, object]:
    activate_runtime_context(
        "encode_regulon_contrast",
        getattr(args, "provenance_overlay_json", None),
    )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    low = float(getattr(args, "low_prevalence", 0.25))
    high = float(getattr(args, "high_prevalence", 0.75))
    assay = getattr(args, "assay", "TF ChIP-seq")
    lib_prefix = getattr(args, "lib_prefix", None) or f"ENCODE_{_safe_name(assay)}_regulon_bgcorrected"

    regulon_zip = getattr(args, "regulon_zip", None)
    regulon_dir = getattr(args, "regulon_dir", None)

    regulon_src_dir: str | None = None
    if regulon_zip:
        src_path = Path(regulon_zip)
        if not src_path.is_file():
            raise FileNotFoundError(f"regulon_zip not found: {src_path}")
        recs = _load_records_from_zip(src_path)
        input_file_rec = input_file_record(str(src_path), "regulon_zip")
    elif regulon_dir:
        src_path = Path(regulon_dir)
        if not src_path.is_dir():
            raise NotADirectoryError(f"regulon_dir not found: {src_path}")
        recs = _load_records_from_dir(src_path)
        input_file_rec = None
        regulon_src_dir = str(src_path)
    else:
        raise ValueError("Provide --regulon_dir or --regulon_zip")

    if not recs:
        raise ValueError("No regulon records loaded from input")

    n = len(recs)

    # Global binding prevalence count: how many factors target each gene
    count: collections.Counter[str] = collections.Counter()
    for _, genes, _ in recs:
        for g in genes:
            count[g] += 1
    universe = set(count)

    n_up = n_down = 0
    up_counts: list[int] = []
    down_counts: list[int] = []

    for setdir_name, genes, target in recs:
        safe_tgt = _safe_name(target)

        # Up: bound by this factor AND LOO-prevalence < low (specifically-bound)
        up: set[str] = set()
        for g in genes:
            loo_prev = (count[g] - 1) / (n - 1) if n > 1 else 0.0
            if loo_prev < low:
                up.add(g)

        # Down: NOT bound here AND LOO-prevalence > high (promiscuous elsewhere)
        down: set[str] = set()
        for g in universe:
            if g in genes:
                continue
            loo_prev = count[g] / (n - 1) if n > 1 else 0.0
            if loo_prev > high:
                down.add(g)

        if up:
            up_name = f"{lib_prefix}_{safe_tgt}_specific_Up"
            _write_contrast_set(
                set_dir=out_dir / up_name,
                set_name=up_name,
                description=(
                    f"Genes specifically targeted by {target} ({assay}) "
                    f"vs cross-factor binding prevalence (LOO, {n} factors, excl. self, "
                    f"< {low}). Removes promiscuous/HOT-region binding. "
                    f"GRCh38; ENCODE/NHGRI; public. {ENCODE_CITATION}"
                ),
                genes=sorted(up),
                assay=assay,
                target=target,
                direction="Up",
                low_threshold=low,
                high_threshold=high,
                n_factors=n,
                input_file_rec=input_file_rec,
                regulon_src_dir=regulon_src_dir,
            )
            up_counts.append(len(up))
            n_up += 1

        if down:
            down_name = f"{lib_prefix}_{safe_tgt}_promiscuous_Down"
            _write_contrast_set(
                set_dir=out_dir / down_name,
                set_name=down_name,
                description=(
                    f"Genes NOT targeted by {target} ({assay}) but commonly bound elsewhere "
                    f"(LOO-prevalence > {high} across {n} factors, excl. self). "
                    f"Broadly bound by other factors but specifically absent here. "
                    f"GRCh38; ENCODE/NHGRI; public. {ENCODE_CITATION}"
                ),
                genes=sorted(down),
                assay=assay,
                target=target,
                direction="Down",
                low_threshold=low,
                high_threshold=high,
                n_factors=n,
                input_file_rec=input_file_rec,
                regulon_src_dir=regulon_src_dir,
            )
            down_counts.append(len(down))
            n_down += 1

    total_counts = up_counts + down_counts
    return {
        "n_factors": n,
        "n_up_sets": n_up,
        "n_down_sets": n_down,
        "n_genes_min": min(total_counts) if total_counts else 0,
        "n_genes_max": max(total_counts) if total_counts else 0,
        "out_dir": str(out_dir),
    }
