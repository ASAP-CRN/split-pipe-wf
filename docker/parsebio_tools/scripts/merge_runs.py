#!/usr/bin/env python3

"""
Merge the cell-typed AnnData objects of each Parse run (kit) into one object.

Runs share the reference annotation, so their genes must be identical; cell IDs already
carry the run (added by merge_and_plot_qc), so they must be unique across runs. X stays
raw counts; uns (per-sample scrublet parameters) is dropped.

Usage
    merge_runs.py --adata-inputs <run1.mmc.h5ad> <run2.mmc.h5ad> ... --adata-output <merged.h5ad>
"""

import argparse
import logging
import anndata as ad
from pathlib import Path

logging.basicConfig(level=logging.INFO)

REQUIRED_OBS = ["batch", "batch_id", "sample", "subject_id", "parse_sample", "cell_type"]


def check_runs(adatas: list[ad.AnnData], paths: list[Path]) -> None:
    """Raise unless the runs can be stacked: same genes, expected obs, distinct batches."""
    for adata, path in zip(adatas, paths):
        if adata.n_obs == 0:
            raise ValueError(f"{path}: no cells")
        if not adata.var_names.equals(adatas[0].var_names):
            raise ValueError(f"{path}: genes differ from {paths[0]}; runs must share the reference")
        missing = [c for c in REQUIRED_OBS if c not in adata.obs]
        if missing:
            raise ValueError(f"{path}: missing obs columns {missing}")

    batches = [frozenset(adata.obs["batch"].unique()) for adata in adatas]
    if len(set(batches)) != len(batches):
        raise ValueError(f"Runs share a batch value: {[sorted(b) for b in batches]}")


def main(args: argparse.Namespace):
    adatas = [ad.read_h5ad(path) for path in args.adata_inputs]
    check_runs(adatas, args.adata_inputs)

    merged = ad.concat(adatas, join="inner", merge="same")
    if not merged.obs_names.is_unique:
        raise ValueError("Cell IDs are not unique across runs")

    for adata, path in zip(adatas, args.adata_inputs):
        logging.info("%s: %d cells", path.name, adata.n_obs)
    logging.info("Merged: %d cells x %d genes", merged.n_obs, merged.n_vars)

    merged.write_h5ad(filename=args.adata_output, compression="lzf")
    logging.info("Saved merged AnnData object to %s", args.adata_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge per-run cell-typed Parse AnnData objects into one object."
    )
    parser.add_argument(
        "--adata-inputs",
        type=Path,
        nargs="+",
        required=True,
        help="Cell-typed AnnData objects, one per Parse run"
    )
    parser.add_argument(
        "--adata-output",
        type=Path,
        required=True,
        help="Output h5ad"
    )

    args = parser.parse_args()
    main(args)
