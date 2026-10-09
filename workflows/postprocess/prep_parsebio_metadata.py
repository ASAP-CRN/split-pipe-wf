#!/usr/bin/env python3
# Direct reference: https://www.sc-best-practices.org/preprocessing-visualization/quality-control/

import argparse
import sys
import logging
import scanpy as sc
import pandas as pd
from anndata import AnnData
from pathlib import Path

logging.basicConfig(level=logging.INFO)


def create_adata_object(data_path: Path) -> AnnData:
    """
    Create an AnnData object from the given data and optional metadata.

    TODO: doc
    Parameters:
    - data: The main data matrix (e.g., gene expression matrix).
    - obs: Optional observation metadata (rows).
    - var: Optional variable metadata (columns).
    - uns: Optional unstructured metadata.

    Returns:
    - AnnData object.
    """
    mtx_path = data_path / "count_matrix.mtx"
    gene_path = data_path / "all_genes.csv"
    cell_meta_path = data_path / "cell_metadata.csv"
    
    if not mtx_path.exists():
        raise FileNotFoundError(f"Matrix file not found: {mtx_path}")
    if not gene_path.exists():
        raise FileNotFoundError(f"Gene file not found: {gene_path}")
    if not cell_meta_path.exists():
        raise FileNotFoundError(f"Cell metadata file not found: {cell_meta_path}")

    adata = sc.read_mtx(mtx_path)
    gene_data = pd.read_csv(gene_path)
    cell_meta_df = pd.read_csv(cell_meta_path)
    
    adata.var_names = gene_data["gene_name"]
    
    # TODO: Here or in QC call
    adata.var_names_make_unique() 
    
    adata.obs = cell_meta_df
    adata.obs_names = adata.obs["bc_wells"]

    logging.info("Created AnnData object with shape %s", adata.shape)
    return adata


def add_vireo_joined_metadata(adata: AnnData, vireo_joined_metadata_path: Path) -> AnnData:
    """
    Add Vireo metadata to the AnnData object.

    TODO: doc
    Parameters:
    - adata: The AnnData object to update.
    - vireo_output_path: Path to the Vireo output directory containing metadata files.

    Returns:
    - Updated AnnData object with Vireo metadata added.
    """
    if not vireo_joined_metadata_path.exists():
        raise FileNotFoundError(f"Vireo metadata file not found: {vireo_joined_metadata_path}")

    vireo_meta_df = pd.read_csv(vireo_joined_metadata_path)
    
    assert all(adata.obs_names == vireo_meta_df["bc_wells"]), \
        "Observation names and Vireo metadata are not aligned"
        
    adata.obs = vireo_meta_df    

    return adata


def add_qc_metrics(adata: AnnData) -> AnnData:
    """
    Add QC metrics to the AnnData object.

    TODO doc
    Parameters:
    - adata: The AnnData object to update.

    Returns:
    - Updated AnnData object with QC metrics added.
    """
    adata.var["mt"] = adata.var_names.str.contains(r"^MT-", case=False, regex=True)
    if not adata.var["mt"].any():
        logging.warning("No mitochondrial genes found using case-insensitive MT- prefix")

    adata.var["rb"] = adata.var_names.str.contains(r"^(RPL|RPS|MRPL|MRPS)", case=False, regex=True)
    if not adata.var["ribo"].any():
        logging.warning("No ribosomal genes found using case-insensitive RPS/RPL/MRPL/MRPS prefix")

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["rb", "mt"],
        percent_top=None,
        log1p=False,  # TODO: check, doc is true, scrnaseq wf is false
        inplace=True,
    )
    
    return adata



def main():
    # TODO: doc 
    parser = argparse.ArgumentParser(description="TODO")
    parser.add_argument("--data_dir", type=Path, help="Path to the directory containing the split-pipe output")    
    parser.add_argument("--vireo_metadata_path", type=Path, help="Path to the Vireo joined metadata CSV file")    
    parser.add_argument("--output", type=Path, help="Path to the output object")
    args = parser.parse_args()

    adata = create_adata_object(data_path=args.data_dir)
    adata = add_vireo_joined_metadata(adata, args.vireo_metadata_path)
    adata.write_h5ad(filename=args.output, compression="gzip")
    logging.info("Saved AnnData object to %s", args.output)


if __name__ == "__main__":
    main()