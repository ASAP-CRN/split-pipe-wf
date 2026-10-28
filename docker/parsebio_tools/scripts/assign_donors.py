#!/usr/bin/env python3

"""
Assign a subject to every cell in one Parse sample's split-pipe cell metadata, from
either a donor demultiplexing table (vireo) or a single known subject (un-pooled sample).
Fails if too few cells match the vireo table, which indicates mismatched inputs (e.g.
sublibrary order differing from the run vireo was computed on).

Output: cell_metadata.csv with split-pipe's "sample" renamed to parse_sample, plus
subject_id ("unassigned" for cells vireo did not call) and vireo_* statistics.

Usage
    python assign_donors.py --cell-meta <cell_metadata.csv> --vireo <vireo.tsv> --output <out.csv>
    python assign_donors.py --cell-meta <cell_metadata.csv> --subject-id <subject_id> --output <out.csv>
"""

import argparse
import logging
import numpy as np
import pandas as pd
from pathlib import Path

logging.basicConfig(level=logging.INFO)


def load_table(file_path: Path, check_join_key: str, sep: str = ",") -> pd.DataFrame:
    """
    Loader for both cell metadata and vireo donor mapping tables that checks that
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


# Team Voet PB-002 pool_5: Contributors confirmed (2026-10-26) that donor9 in the
# Vireo output is ASA_137. Since Vireo generically labels donorN, this replacement
# is keyed to be specific to the dataset and Parse sample ("pool_5"). This name
# exists only in PB-002 and does not/will never exist in PB-001.
VIREO_DONOR_RENAMES = {
    ("team-voet-pmdbs-sn-rnaseq-parsebio", "pool_5"): {"donor9": "ASA_137"},
}


def rename_vireo_donors(join_df: pd.DataFrame, dataset_id: str) -> pd.DataFrame:
    """
    Apply confirmed vireo donor -> subject renames for this dataset and Parse sample
    (VIREO_DONOR_RENAMES). Raises if an expected vireo donor is absent or the target subject
    already has vireo-assigned cells, i.e. if the inputs no longer match the confirmation.
    """
    parse_sample = join_df["parse_sample"].iloc[0]
    renames = VIREO_DONOR_RENAMES.get((dataset_id, parse_sample), {})
    for vireo_donor, subject_id in renames.items():
        is_donor = join_df["subject_id"] == vireo_donor
        if not is_donor.any():
            raise ValueError(f"Expected vireo donor {vireo_donor} in {dataset_id} {parse_sample}; not found")
        if (join_df["subject_id"] == subject_id).any():
            raise ValueError(f"{subject_id} already has vireo-assigned cells in {parse_sample}; not renaming {vireo_donor}")
        logging.info(
            "Renaming vireo donor %s to %s for %d cells (%s %s)",
            vireo_donor, subject_id, is_donor.sum(), dataset_id, parse_sample,
        )
        join_df.loc[is_donor, "subject_id"] = subject_id
        
    # Any vireo donor still unnamed has no subject to join
    leftover = sorted(join_df.loc[join_df["subject_id"].str.fullmatch(r"donor\d+", na=False), "subject_id"].unique())
    if leftover:
        raise ValueError(
            f"Vireo donors without a subject in {dataset_id} {parse_sample}: {leftover}"
        )

    return join_df


def assign_subjects(
    cell_meta: pd.DataFrame,
    vireo_df: pd.DataFrame | None,
    subject_id: str | None,
    dataset_id: str,
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
        join_df = join_df.rename(columns={"donor_id": "subject_id", **vireo_renames})

        # Confirmed vireo donorN -> subject_id rename
        join_df = rename_vireo_donors(join_df, dataset_id)
        
        # Check the proportion of missing subject IDs and set NAs to "unassigned"
        check_missing_subject_ids(join_df)
        join_df["subject_id"] = join_df["subject_id"].fillna("unassigned")
        
    else:
        join_df = cell_meta.copy()
        join_df["subject_id"] = subject_id
        for col in vireo_renames.values():
            join_df[col] = np.nan

    return join_df


def check_missing_subject_ids(
    join_df: pd.DataFrame, min_prop: float = 0.05, warn_prop: float = 0.5, min_sublib_prop: float = 0.05
) -> None:
    """
    Check the proportion of cells with a vireo subject. Below min_prop the cell metadata
    and vireo table do not describe the same cells (a random overlap is ~1-3%, e.g. the
    wrong vireo file) and this raises. Any single sublibrary below min_sublib_prop also
    raises: a sublibrary order differing from the vireo run leaves the overall proportion
    high but drops the swapped sublibraries to random overlap. Below warn_prop it warns:
    split-pipe's cell calls can legitimately exceed those vireo was run on (PB-002 pool_1,
    where the contributors' run called far fewer cells in most sublibraries; lowest
    sublibrary observed 0.152).
    """
    non_na_proportion = join_df["subject_id"].notna().mean()
    by_sublib = (
        join_df.assign(sublib=join_df["bc_wells"].str.split("__").str[1])
        .groupby("sublib")["subject_id"]
        .apply(lambda s: round(s.notna().mean(), 3))
        .to_dict()
    )
    low_sublibs = {sublib: prop for sublib, prop in by_sublib.items() if prop < min_sublib_prop}
    if non_na_proportion < min_prop:
        raise ValueError(
            f"Proportion of cells with a vireo subject is {non_na_proportion:.3f}, "
            f"below {min_prop}; per sublibrary: {by_sublib}"
        )
    elif low_sublibs:
        raise ValueError(
            f"Sublibraries with a vireo subject proportion below {min_sublib_prop} "
            f"(sublibrary order may differ from the vireo run): {low_sublibs}; per sublibrary: {by_sublib}"
        )
    elif non_na_proportion < warn_prop:
        logging.warning(
            "Proportion of cells with a vireo subject is low: %.3f; per sublibrary: %s",
            non_na_proportion, by_sublib,
        )
    else:
        logging.info(
            "Proportion of cells with a vireo subject: %.3f; per sublibrary: %s",
            non_na_proportion, by_sublib,
        )


def main(args: argparse.Namespace):

    # split-pipe meta output as .csv, provided vireo tables are .tsv
    cell_meta = load_table(args.cell_meta, sep=",", check_join_key="bc_wells")

    if args.vireo is not None:
        vireo_df = load_table(args.vireo, sep="\t", check_join_key="cell")
    else:
        vireo_df = None

    join_data = assign_subjects(
        cell_meta,
        vireo_df,
        subject_id=args.subject_id,
        dataset_id=args.dataset_id
    )
    
    join_data.to_csv(args.output, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Assign donors to split-pipe cell metadata from a vireo table or a single known subject."
    )
    parser.add_argument(
        "--cell-meta",
        type=Path,
        required=True,
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
    parser.add_argument(
        "--dataset-id",
        required=True,
        help="ASAP dataset ID; scopes confirmed vireo donor renames (VIREO_DONOR_RENAMES)"
    )

    args = parser.parse_args()
    main(args)
