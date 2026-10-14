#!/usr/bin/env python3

"""
Join the donor/subject IDs and additional vireo statistics to sample's cell metadata 
table. Warn if a significant proportion of donor IDs are missing, which may indicate 
a mismatch between the provided tables.

Usage
    python join_donor_ids_from_vireo.py 
        --cell_meta <cell_meta.csv> \
        --vireo <vireo.tsv> \
        --subject_id <subject_id> \
        --output <output.csv>
"""

import argparse
import pandas as pd
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)


def load_table(file_path: Path, check_join_key: str, sep: str = ",") -> pd.DataFrame:
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
    cell_meta: pd.DataFrame, 
    vireo_df: pd.DataFrame | None, 
    subject_id: str | None,
    left_on: str = "bc_wells",
    right_on: str = "cell",
    ) -> pd.DataFrame:
    """
    If vireo_df is provided, join to the cell metadata table to add donor/subject
    IDs and additional vireo statistics. The Vireo file calls it the donor_id, 
    rename to subject_id to match ASAP metadata. If subject_id is provided, add 
    this information directly and fill out the remaining vireo columns as NA.
    
    This branching is because for the Team Voet Parse Bio dataset, every subject
    but one (ASA_063) was genetically pooled and required vireo-based donor assignment.
    For the one subject we can just plug it in directly while the remaining vireo 
    columns will be filled with NA.
    """
    if vireo_df is not None:
        join_df = pd.merge(
            cell_meta,
            vireo_df,
            left_on=left_on,
            right_on=right_on,
            how="left"
        )
        join_df.drop(columns=["cell"], inplace=True)
        join_df.rename(columns={"donor_id": "subject_id"}, inplace=True)
        
        # Check the proportion of missing subject IDs and set NAs to "unassigned"
        check_missing_subject_ids(join_df)
        join_df["subject_id"] = join_df["subject_id"].fillna("unassigned")
        
    else:
        join_df = cell_meta.copy()
        join_df["subject_id"] = subject_id
        vireo_cols = ["prob_max", "prob_doublet", "n_vars", "best_singlet", "best_doublet", "doublet_logLikRatio"]
        for col in vireo_cols:
            join_df[col] = pd.NA
            
    return join_df


def check_missing_subject_ids(join_df: pd.DataFrame, min_prop: float = 0.8) -> None:
    """
    Check the proportion of missing subject IDs and log a warning if it exceeds the threshold.
    """
    non_na_proportion = join_df['subject_id'].notna().mean()
    if non_na_proportion < min_prop:
        logging.error(f"Proportion of non-NA subject_id is low: {non_na_proportion}")
        raise ValueError(f"Proportion of non-NA subject_id is below the minimum threshold of {min_prop}: {non_na_proportion}")
    else:
        logging.info(f"Proportion of non-NA subject_id: {non_na_proportion}")


def main():
    parser = argparse.ArgumentParser(description="Join donor/subject IDs from Vireo to cell metadata.")
    parser.add_argument("--cell_meta", type=Path, help="Path to the cell metadata CSV file.")
    
    # Either a Vireo file or a subject ID must be provided, but not both or neither
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--vireo", type=Path, help="Path to the Vireo TSV file.")
    group.add_argument("--subject_id", type=str, help="Donor/subject ID for the sample.")
    
    parser.add_argument("--output", type=Path, help="Path to the output CSV file.")
    args = parser.parse_args()

    # split-pipe meta output as .csv, provided vireo tables are .tsv
    cell_meta = load_table(args.cell_meta, sep=",", check_join_key="bc_wells")
    
    if args.vireo is not None:
        vireo_df = load_table(args.vireo, sep="\t", check_join_key="cell")
    else:
        vireo_df = None

    join_data = join_cell_meta_and_vireo(cell_meta, vireo_df, subject_id=args.subject_id)

    join_data.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()
