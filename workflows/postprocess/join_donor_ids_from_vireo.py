#!/usr/bin/env python3

"""
Join the donor/subject IDs and additional vireo statistics from the provided Team
Voet Vireo table to the parse sample's cell metadata table to demultiplex cells 
by donor. Warn if a significant proportion of donor IDs are missing, which may indicate 
a mismatch between the provided tables.

Usage
    python join_donor_ids_from_vireo.py 
        --cell-meta <cell_meta.csv> \
        --vireo <vireo.tsv> \
        --subject-id <subject_id> \
        --output <output.csv>
"""

import argparse
import logging
import numpy as np
import pandas as pd
from pathlib import Path

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


# Vireo per-cell statistics: prefix so they are not confused with scrublet's doublet_score
# or AnnData's n_vars (vireo's n_vars is the number of informative SNPs)
VIREO_STAT_COLS = ["prob_max", "prob_doublet", "n_vars", "best_singlet", "best_doublet", "doublet_logLikRatio"]


def join_cell_meta_and_vireo(
    cell_meta: pd.DataFrame, 
    vireo_df: pd.DataFrame | None, 
    subject_id: str | None,
) -> pd.DataFrame:
    """
    Assign a subject to every cell in one Parse sample's cell metadata.

    - vireo_df (pooled sample): left-join vireo's donor calls; vireo's donor_id becomes
      subject_id, and cells vireo did not call are set to "unassigned".
    - subject_id (un-pooled sample): every cell gets that subject; vireo columns are NA.

    split-pipe's "sample" column is the Parse sample (e.g. pool_1), renamed to parse_sample
    so "sample" is free for the ASAP sample ID. Vireo statistics are prefixed with "vireo_".
    """
    cell_meta = cell_meta.rename(columns={"sample": "parse_sample"})
    vireo_renames = {c: f"vireo_{c}" for c in VIREO_STAT_COLS}
    
    if vireo_df is not None:
        join_df = pd.merge(cell_meta, vireo_df, left_on="bc_wells", right_on="cell", how="left")
        join_df.drop(columns=["cell"], inplace=True)
        join_df.rename(columns={"donor_id": "subject_id", **vireo_renames})
        
        # Check the proportion of missing subject IDs and set NAs to "unassigned"
        check_missing_subject_ids(join_df)
        join_df["subject_id"] = join_df["subject_id"].fillna("unassigned")
        
    else:
        join_df = cell_meta.copy()
        join_df["subject_id"] = subject_id
        for col in vireo_renames.values():
            join_df[col] = np.nan
            
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


def main(args: argparse.Namespace):

    # split-pipe meta output as .csv, provided vireo tables are .tsv
    cell_meta = load_table(args.cell_meta, sep=",", check_join_key="bc_wells")
    
    if args.vireo is not None:
        vireo_df = load_table(args.vireo, sep="\t", check_join_key="cell")
    else:
        vireo_df = None

    join_data = join_cell_meta_and_vireo(cell_meta, vireo_df, subject_id=args.subject_id)
    join_data.to_csv(args.output, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Assign donors to split-pipe cell metadata from a vireo table or a single known subject."
    )
    parser.add_argument(
        "--cell-meta",
        type=Path,
        help="split-pipe cell_metadata.csv"
    )
    # Either a vireo table (pooled Parse sample) or a single subject ID (un-pooled), not both
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--vireo",
        type=Path,
        help="Vireo donor assignment TSV"
    )
    group.add_argument(
        "--subject-id",
        type=str,
        help="Subject ID for an un-pooled Parse sample"
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output CSV"
    )

    args = parser.parse_args()
    main(args)
