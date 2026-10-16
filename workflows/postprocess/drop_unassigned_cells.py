#!/usr/bin/env python3

"""
Drop cells without an ASAP sample (vireo doublet/unassigned) from the merged AnnData
object, before sc-rnaseq-wf's filter. Removed cells are written in the same layout as
filter's --output-removed-cells-file so the two reports can be stacked.

Usage
    python drop_unassigned_cells.py \
        --adata-input <merged.h5ad> \
        --adata-output <assigned.h5ad> \
        --output-removed-cells-file <unassigned_removed_cells.csv>
"""

import argparse
import logging
from pathlib import Path

import pandas as pd
import scanpy as sc
from anndata import AnnData

logging.basicConfig(level=logging.INFO)

# Same metadata columns filter carries into its removed-cells report
REPORT_COLS = ["sample", "batch_id", "team", "dataset", "brain_region_level_1",
               "pct_counts_mt", "doublet_score", "total_counts", "n_genes_by_counts"]


def drop_unassigned(adata_input: Path) -> tuple[AnnData, pd.DataFrame]:
    """
    Keep cells that resolved to an ASAP sample. Cells without one are vireo doublet or
    unassigned calls, including cells split-pipe called that vireo did not.
    """
    adata_backed = sc.read_h5ad(adata_input, backed="r")
    obs = adata_backed.obs
    keep = obs["sample"].notna()

    calls = obs.loc[~keep, "subject_id"].value_counts()
    logging.info("Keeping %d/%d cells; removed by vireo call: %s",
                 keep.sum(), len(keep), calls.to_dict())

    cols = [c for c in REPORT_COLS + ["parse_sample", "subject_id"] if c in obs.columns]
    removed_cells = obs.loc[~keep, cols].copy()
    removed_cells["failed_vireo_assignment"] = True

    adata = adata_backed[keep.to_numpy()].to_memory()
    adata_backed.file.close()

    return adata, removed_cells


def main(args: argparse.Namespace):
    adata, removed_cells = drop_unassigned(args.adata_input)
    adata.write_h5ad(filename=args.adata_output, compression="lzf")
    removed_cells.to_csv(args.output_removed_cells_file, index=True, index_label="barcode")
    logging.info("Wrote %s and %s", args.adata_output, args.output_removed_cells_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Drop vireo doublet/unassigned cells before filter")
    parser.add_argument(
        "--adata-input",
        type=Path,
        required=True,
        help="Merged AnnData object from merge_and_plot_qc"
    )
    parser.add_argument(
        "--adata-output",
        type=Path,
        required=True,
        help="AnnData object with vireo non-singlets removed"
    )
    parser.add_argument(
        "--output-removed-cells-file",
        type=Path,
        required=True,
        help="CSV of removed cells, same layout as filter's removed-cells report"
    )
    args = parser.parse_args()
    main(args)
