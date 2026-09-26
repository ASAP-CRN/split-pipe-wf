#!/usr/bin/env python3

"""
Join the donor IDs and additional vireo statistics to sample's cell metadata table.
Warn if a significant proportion of donor IDs are missing, which may indicate a
mismatch between the provided tables.

Usage
    python join_donor_ids_from_vireo.py <cell_meta.csv> <vireo.tsv> <output.csv>
"""

import argparse
import pandas as pd
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)


def load_table(file_path: Path, check_join_key: str, sep: str = ","):
    """
    Loader for both cell metadata and vireo donoor mapping tables that checks that
    the expected join key is present in the table.
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"File {file_path} does not exist")
    
    df = pd.read_csv(file_path, sep=sep)
    
    if check_join_key not in df.columns:
        raise ValueError(f"Join key {check_join_key} not found in {file_path}")
    
    return df


def join_cell_meta_and_vireo(
    sample_meta: pd.DataFrame, 
    vireo_df: pd.DataFrame, 
    left_on: str = "bc_wells", 
    right_on: str = "cell"
    ) -> pd.DataFrame:
    """
    Join adds the donor information and additional statistics from the vireo mapping
    to the cell metadata.
    """
    join_df = pd.merge(
        sample_meta,
        vireo_df,
        left_on=left_on,
        right_on=right_on,
        how="left"
    )
    join_df.drop(columns=["cell"], inplace=True)
    
    return join_df


def check_missing_donor_ids(join_df: pd.DataFrame, min_prop: float = 0.5) -> None:
    """
    Check the proportion of missing donor IDs and log a warning if it exceeds the threshold.
    """
    non_na_proportion = join_df['donor_id'].notna().mean()
    if non_na_proportion < min_prop:
        logging.warning(f"Proportion of non-NA donor_id is low: {non_na_proportion}")
    else:
        logging.info(f"Proportion of non-NA donor_id: {non_na_proportion}")


def main():
    parser = argparse.ArgumentParser(description="Join donor IDs from Vireo to cell metadata.")
    parser.add_argument("--cell_meta", type=Path, help="Path to the cell metadata CSV file.")
    parser.add_argument("--vireo", type=Path, help="Path to the Vireo TSV file.")
    parser.add_argument("--output", type=Path, help="Path to the output CSV file.")
    args = parser.parse_args()

    # split-pipe meta output as .csv, provided vireo tables are .tsv
    cell_meta = load_table(args.cell_meta, sep=",", check_join_key="bc_wells")
    vireo_df = load_table(args.vireo, sep="\t", check_join_key="cell")

    join_data = join_cell_meta_and_vireo(cell_meta, vireo_df)
    check_missing_donor_ids(join_data)

    join_data.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()
