#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess


downloaded_files = "/mnt/var/single_cell"
s3_in = os.environ["INPUT_PATH"]
s3_out = os.environ["OUTPUT_PATH"]


def download(dataset: str, cell_type: str) -> None:
    path = f"{s3_in}/out/single_cell/staging/h5ad/{dataset}/{cell_type}.h5ad"
    subprocess.check_call(["aws", "s3", "cp", path, "inputs/data.h5ad"])


def run_liger() -> None:
    cmd = [
        "Rscript",
        f"{downloaded_files}/inmf_liger_v2-4_mod.R",
        "inputs/data.h5ad",
        "outputs",
        "donor_id",
        "cell_type__kp",
        "50000",
        "50",
    ]
    subprocess.check_call(cmd)


def upload(dataset: str) -> None:
    if os.path.exists("./outputs/"):
        path = f"{s3_out}/out/single_cell/staging/liger/{dataset}/"
        subprocess.check_call(["aws", "s3", "cp", "./outputs/", path, "--recursive"])
        shutil.rmtree("outputs")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True, type=str, help="Dataset name")
    parser.add_argument("--cell-type", required=True, type=str, help="Cell type name")
    args = parser.parse_args()

    download(args.dataset, args.cell_type)
    run_liger()
    upload(args.dataset)
    shutil.rmtree("inputs")


if __name__ == "__main__":
    main()
