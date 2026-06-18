#!/usr/bin/env python3
import argparse
import gzip
import math
import os
import shutil
import subprocess

import anndata as ad
import numpy as np
from scipy.sparse import csc_matrix, vstack


s3_in = os.environ["INPUT_PATH"]
s3_out = os.environ["OUTPUT_PATH"]


def download_data(dataset: str) -> None:
    subprocess.check_call(["aws", "s3", "cp", f"{s3_in}/single_cell/{dataset}/norm_counts.tsv.gz", "inputs/"])
    subprocess.check_call(["aws", "s3", "cp", f"{s3_in}/single_cell/{dataset}/sample_metadata.tsv.gz", "inputs/"])


def get_metadata_maps():
    cell_type_map = {}
    donor_map = {}
    ncount_map = {}
    with gzip.open("inputs/sample_metadata.tsv.gz", "rt") as f:
        header = f.readline().strip().split("\t")
        for line in f:
            json_line = dict(zip(header, line.strip().split("\t")))
            n_count = json_line["QC:nCount_RNA"]
            if len(n_count) > 0 and float(n_count) % 1 == 0:
                cell_type_map[json_line["ID"]] = json_line["cell_type__kp"].replace("/", " ")
                donor_map[json_line["ID"]] = json_line["DI:Dataset"]
                ncount_map[json_line["ID"]] = int(n_count)
    return cell_type_map, donor_map, ncount_map


def get_sparse_array(cell_type_map, donor_map, ncount_map):
    cell_types = list(set(cell_type_map.values()))
    with gzip.open("inputs/norm_counts.tsv.gz", "rt") as f:
        cells = f.readline().strip().split("\t")[1:]
        p_map = {
            cell_type: 100000 / len([cell for cell in cells if cell_type_map.get(cell, "") == cell_type])
            for cell_type in cell_types
        }
        cells_idx_dict = {
            cell_type: [
                (idx, ncount_map[cell])
                for idx, cell in enumerate(cells)
                if cell_type_map.get(cell, "") == cell_type and np.random.rand() < p_map[cell_type]
            ]
            for cell_type in cell_types
        }
        gene, data = f.readline().strip().split("\t", 1)
        genes = [gene]
        formatted_data = list(map(float, data.split("\t")))
        A_dict = {
            cell_type: csc_matrix(
                [[int(round((math.exp(formatted_data[idx]) - 1) * ncount / 1e4)) for idx, ncount in cells_idx_dict[cell_type]]]
            )
            for cell_type in cell_types
        }
        count = 1
        B_dict = {cell_type: [] for cell_type in cell_types}
        for line in f:
            if count == 100:
                for cell_type in cell_types:
                    A_dict[cell_type] = vstack([A_dict[cell_type], csc_matrix(B_dict[cell_type])])
                B_dict = {cell_type: [] for cell_type in cell_types}
                count = 0
            gene, data = line.strip().split("\t", 1)
            if gene not in genes:
                line_to_append = list(map(float, data.split("\t")))
                for cell_type in cell_types:
                    B_dict[cell_type].append(
                        [int(round((math.exp(line_to_append[idx]) - 1) * ncount / 1e4)) for idx, ncount in cells_idx_dict[cell_type]]
                    )
                count += 1
                genes.append(gene)
        for cell_type in cell_types:
            A_dict[cell_type] = vstack([A_dict[cell_type], csc_matrix(B_dict[cell_type])])
    return {
        cell_type: ad.AnnData(
            A_dict[cell_type].T,
            obs={
                "obs_names": [cells[idx] for idx, _ in cells_idx_dict[cell_type]],
                "cell_type__kp": [cell_type_map[cells[idx]] for idx, _ in cells_idx_dict[cell_type]],
                "donor_id": [donor_map[cells[idx]] for idx, _ in cells_idx_dict[cell_type]],
            },
            var={"var_names": genes},
        )
        for cell_type in cell_types
    }


def upload(dataset: str, cell_type_adata_map) -> None:
    for cell_type, adata in cell_type_adata_map.items():
        adata.write_h5ad(f"{cell_type}.h5ad")
        subprocess.check_call(["aws", "s3", "cp", f"{cell_type}.h5ad", f"{s3_out}/out/single_cell/staging/h5ad/{dataset}/"])
        os.remove(f"{cell_type}.h5ad")
    shutil.rmtree("inputs")


def run() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True, type=str, help="Dataset name")
    args = parser.parse_args()

    download_data(args.dataset)
    cell_type_map, donor_map, ncount_map = get_metadata_maps()
    cell_type_adata_dict = get_sparse_array(cell_type_map, donor_map, ncount_map)
    upload(args.dataset, cell_type_adata_dict)


if __name__ == "__main__":
    run()
