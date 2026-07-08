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

_Record = tuple[str, frozenset[str], str, str]
# (setdir_label, genes, biosample_label, group)


def _safe_name(s: str, max_len: int = 60) -> str:
    return "".join(c if c.isalnum() else "_" for c in s)[:max_len]


def _load_valid_symbols(tsv_path: Path) -> set[str] | None:
    """Load valid HGNC symbols from the first column of a TSV (skip header)."""
    symbols: set[str] = set()
    with tsv_path.open("r", encoding="utf-8") as fh:
        next(fh, None)  # skip header
        for line in fh:
            sym = line.split("\t")[0].strip()
            if sym:
                symbols.add(sym)
    return symbols if symbols else None


def _load_records_from_dir(
    src: Path,
    group_key: str | None,
    lib_filter: str | None,
    valid: set[str] | None,
) -> list[_Record]:
    """Load (label, genes, biosample, group) from a directory of accessible-gene set subdirs."""
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
        if valid is not None:
            genes = {g for g in genes if g in valid}
        if not genes:
            continue
        label = sub.name
        group = "ALL"
        meta_path = sub / "geneset.meta.json"
        if meta_path.is_file():
            try:
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
                label = meta.get("biosample", label) or label
                if lib_filter and lib_filter not in str(meta.get("library", "")):
                    continue
                if group_key:
                    group = str(meta.get(group_key, "ALL"))
            except Exception:
                if lib_filter:
                    continue
        recs.append((sub.name, frozenset(genes), label, group))
    return recs


def _load_records_from_zip(
    zpath: Path,
    group_key: str | None,
    lib_filter: str | None,
    valid: set[str] | None,
) -> list[_Record]:
    """Load records from a zip file of per-biosample accessible-gene set dirs."""
    recs: list[_Record] = []
    with zipfile.ZipFile(zpath, "r") as zf:
        names = zf.namelist()
        for n in names:
            if not n.endswith("geneset.tsv"):
                continue
            setdir = n.split("/")[-2]
            raw = zf.read(n).decode("utf-8", "replace").splitlines()
            genes: set[str] = {ln.strip() for ln in raw[1:] if ln.strip()}
            if valid is not None:
                genes = {g for g in genes if g in valid}
            if not genes:
                continue
            label = setdir
            group = "ALL"
            meta_name = n.rsplit("/", 1)[0] + "/geneset.meta.json"
            if meta_name in names:
                try:
                    meta = json.loads(zf.read(meta_name).decode("utf-8", "replace"))
                    label = meta.get("biosample", label) or label
                    if lib_filter and lib_filter not in str(meta.get("library", "")):
                        continue
                    if group_key:
                        group = str(meta.get(group_key, "ALL"))
                except Exception:
                    if lib_filter:
                        continue
            recs.append((setdir, frozenset(genes), label, group))
    return recs


def _write_contrast_set(
    set_dir: Path,
    set_name: str,
    description: str,
    genes: list[str],
    assay: str,
    biosample: str,
    group: str,
    direction: str,
    low_threshold: float,
    high_threshold: float,
    n_background: int,
    input_file_rec: dict | None,
    accessible_src_dir: str | None = None,
) -> None:
    set_dir.mkdir(parents=True, exist_ok=True)
    with (set_dir / "geneset.tsv").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("gene\n")
        for g in genes:
            fh.write(f"{g}\n")
    write_gmt([(set_name, genes)], set_dir / "genesets.gmt")

    params: dict = {
        "assay": assay,
        "biosample": biosample,
        "group": group,
        "direction": direction,
        "low_prevalence_threshold": low_threshold,
        "high_prevalence_threshold": high_threshold,
        "background_method": "leave_one_out_prevalence",
        "background_n_experiments": n_background,
        "encode_citation": ENCODE_CITATION,
        "caveat": (
            "Peak-call prevalence background (relative accessibility specificity); "
            "not a GC/library-normalized read-count differential (no DESeq2/edgeR)."
        ),
    }
    if accessible_src_dir is not None:
        params["accessible_src_dir"] = accessible_src_dir
    meta = make_metadata(
        converter_name="encode_accessibility_contrast",
        parameters=params,
        data_type="chromatin_accessibility",
        assay=assay,
        organism="human",
        genome_build="GRCh38",
        files=[input_file_rec] if input_file_rec is not None else [],
        gene_annotation={
            "mode": "inherited",
            "source": "encode_accessible_genes_output",
        },
        weights={
            "weight_type": "nonnegative",
            "normalization": {"method": "none", "target_sum": None},
            "aggregation": "loo_prevalence_threshold",
        },
        summary={
            "n_input_features": n_background,
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
        "encode_accessibility_contrast",
        getattr(args, "provenance_overlay_json", None),
    )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    low = float(getattr(args, "low_prevalence", 0.25))
    high = float(getattr(args, "high_prevalence", 0.75))
    assay = getattr(args, "assay", "ATAC-seq")
    group_key = getattr(args, "group_by_key", None) or None
    lib_filter = getattr(args, "lib_filter", None) or None

    # Optional HGNC symbol universe filter
    sym_tsv = getattr(args, "symbol_universe_tsv", None)
    valid: set[str] | None = None
    if sym_tsv:
        valid = _load_valid_symbols(Path(sym_tsv))

    # Load records from directory or zip
    accessible_zip = getattr(args, "accessible_genes_zip", None)
    accessible_dir = getattr(args, "accessible_genes_dir", None)

    accessible_src_dir: str | None = None
    if accessible_zip:
        src_path = Path(accessible_zip)
        if not src_path.is_file():
            raise FileNotFoundError(f"accessible_genes_zip not found: {src_path}")
        recs = _load_records_from_zip(src_path, group_key, lib_filter, valid)
        input_file_rec = input_file_record(str(src_path), "accessible_genes_zip")
    elif accessible_dir:
        src_path = Path(accessible_dir)
        if not src_path.is_dir():
            raise NotADirectoryError(f"accessible_genes_dir not found: {src_path}")
        recs = _load_records_from_dir(src_path, group_key, lib_filter, valid)
        input_file_rec = None
        accessible_src_dir = str(src_path)
    else:
        raise ValueError("Provide --accessible_genes_dir or --accessible_genes_zip")

    if not recs:
        raise ValueError("No accessible-gene records loaded from input")

    # Group records for per-group background computation
    groups: dict[str, list[_Record]] = collections.defaultdict(list)
    for rec in recs:
        groups[rec[3]].append(rec)

    n_up = n_down = 0
    up_counts: list[int] = []
    down_counts: list[int] = []

    for grp, members in sorted(groups.items()):
        n = len(members)
        # Per-group prevalence count: how many biosamples in this group have each gene
        count: collections.Counter[str] = collections.Counter()
        for _, genes, _, _ in members:
            for g in genes:
                count[g] += 1
        universe = set(count)
        prefix = _safe_name(assay)
        if grp != "ALL":
            prefix = f"{prefix}_{_safe_name(grp)}"

        for setdir_label, genes, biosample, _ in members:
            safe_bs = _safe_name(biosample)

            # Up: accessible here AND LOO-prevalence < low
            up: set[str] = set()
            for g in genes:
                loo_prev = (count[g] - 1) / (n - 1) if n > 1 else 0.0
                if loo_prev < low:
                    up.add(g)

            # Down: absent here AND LOO-prevalence > high
            down: set[str] = set()
            for g in universe:
                if g in genes:
                    continue
                loo_prev = count[g] / (n - 1) if n > 1 else 0.0
                if loo_prev > high:
                    down.add(g)

            grp_label = "" if grp == "ALL" else f" [{grp}]"

            if up:
                up_name = f"ENCODE_{prefix}_{safe_bs}_accessible_Up"
                _write_contrast_set(
                    set_dir=out_dir / up_name,
                    set_name=up_name,
                    description=(
                        f"Genes specifically accessible ({assay}{grp_label}) in '{biosample}' "
                        f"vs LOO background (prevalence across {n} experiments, excl. self, "
                        f"< {low}). "
                        f"Removes constitutively-open housekeeping genes. "
                        f"GRCh38; ENCODE/NHGRI; public. {ENCODE_CITATION}"
                    ),
                    genes=sorted(up),
                    assay=assay,
                    biosample=biosample,
                    group=grp,
                    direction="Up",
                    low_threshold=low,
                    high_threshold=high,
                    n_background=n,
                    input_file_rec=input_file_rec,
                    accessible_src_dir=accessible_src_dir,
                )
                up_counts.append(len(up))
                n_up += 1

            if down:
                down_name = f"ENCODE_{prefix}_{safe_bs}_inaccessible_Down"
                _write_contrast_set(
                    set_dir=out_dir / down_name,
                    set_name=down_name,
                    description=(
                        f"Genes specifically inaccessible ({assay}{grp_label}) in '{biosample}' "
                        f"vs LOO background (prevalence across {n} experiments, excl. self, "
                        f"> {high}). "
                        f"Broadly open elsewhere but specifically closed here. "
                        f"GRCh38; ENCODE/NHGRI; public. {ENCODE_CITATION}"
                    ),
                    genes=sorted(down),
                    assay=assay,
                    biosample=biosample,
                    group=grp,
                    direction="Down",
                    low_threshold=low,
                    high_threshold=high,
                    n_background=n,
                    input_file_rec=input_file_rec,
                    accessible_src_dir=accessible_src_dir,
                )
                down_counts.append(len(down))
                n_down += 1

    total = n_up + n_down
    all_counts = up_counts + down_counts
    return {
        "n_groups": total,
        "n_up_sets": n_up,
        "n_down_sets": n_down,
        "n_genes_min": min(all_counts) if all_counts else 0,
        "n_genes_max": max(all_counts) if all_counts else 0,
        "out_dir": str(out_dir),
    }
