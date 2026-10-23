#!/usr/bin/env python3

"""
Build an AnnData object for one Parse sample from split-pipe DGE output: attach donor
assignments, compute QC metrics and Scrublet scores, and join ASAP SAMPLE/SUBJECT metadata.

Input contract with assign_donors.py: a CSV of split-pipe's cell_metadata.csv plus a
subject_id column (and vireo_* columns, NA when the sample was not pooled), with bc_wells
in the same order as count_matrix.mtx rows.

QC metric definitions and the Scrublet call match sc-rnaseq-wf's prep_metadata:
https://github.com/ASAP-CRN/sc-rnaseq-wf/blob/main/docker/sc_tools/scripts/main/prep_metadata
"""

import argparse
import logging
import scanpy as sc
import pandas as pd
from anndata import AnnData
from pathlib import Path

logging.basicConfig(level=logging.INFO)


def create_adata_object(data_path: Path) -> AnnData:
    """
    Load the split-pipe matrix, gene annotations, and cell metadata into AnnData.

    Expects count_matrix.mtx (cells x genes), all_genes.csv, and cell_metadata.csv in
    data_path. Gene symbols (made unique) are the var index, with Ensembl IDs in
    var["gene_id"]. bc_wells is the obs index.
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

    adata.var = gene_data.set_index("gene_name", drop=False)
    adata.var_names_make_unique()
    adata.var.index.name = "gene_symbol_unique"

    adata.obs = cell_meta_df.set_index("bc_wells", drop=False)

    logging.info("Created AnnData object with shape %s", adata.shape)
    return adata


def add_donor_assignments(adata: AnnData, donor_assignments_path: Path) -> AnnData:
    """
    Replace obs with the donor-assigned cell metadata from assign_donors.py.

    The CSV is split-pipe's cell_metadata.csv plus subject_id and vireo_* columns; its
    bc_wells values and order must match the count matrix rows exactly.
    """
    if not donor_assignments_path.exists():
        raise FileNotFoundError(f"Donor assignments file not found: {donor_assignments_path}")

    assignments = pd.read_csv(donor_assignments_path, index_col="bc_wells")

    if not adata.obs_names.equals(assignments.index):
        raise ValueError("Cell order in the donor assignments does not match count_matrix.mtx rows")

    adata.obs = assignments

    return adata


def add_qc_metrics(adata: AnnData) -> AnnData:
    """
    Calculate gene- and cell-level QC metrics, then add Scrublet doublet scores.

    Mitochondrial and ribosomal genes are identified from gene-symbol prefixes, as in
    sc-rnaseq-wf. split-pipe's gene_count and tscp_count are checked against Scanpy's
    n_genes_by_counts and total_counts, then dropped as redundant. Scrublet runs on all
    called cells, including vireo doublet/unassigned, with default parameters;
    predicted_doublet is dropped in favour of doublet_score.
    """
    adata.var["mt"] = adata.var_names.str.contains(r"^MT-", case=False, regex=True)
    if not adata.var["mt"].any():
        logging.warning("No mitochondrial genes found using case-insensitive MT- prefix")

    adata.var["rb"] = adata.var_names.str.contains(r"^(?:RPL|RPS|MRPL|MRPS)", case=False, regex=True)
    if not adata.var["rb"].any():
        logging.warning("No ribosomal genes found using case-insensitive RPS/RPL/MRPL/MRPS prefix")

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["rb", "mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    # Remove redundant columns with split-pipe output
    if not (adata.obs["gene_count"] == adata.obs["n_genes_by_counts"]).all():
        raise ValueError("split-pipe gene_count disagrees with n_genes_by_counts")
    if not (adata.obs["tscp_count"] == adata.obs["total_counts"]).all():
        raise ValueError("split-pipe tscp_count disagrees with total_counts")
    adata.obs.drop(columns=["gene_count", "tscp_count"], inplace=True)

    # Scrublet: drop predicted_doublet; rely on the probability of doublet (i.e., doublet_score)
    sc.pp.scrublet(adata)
    adata.obs.drop("predicted_doublet", axis=1, inplace=True)

    return adata


def build_sample_lookup(
    sample_df: pd.DataFrame,
    subject_df: pd.DataFrame,
    batch: str,
    sample_map: pd.DataFrame | None = None,  # optional: parse_sample, subject_id, sample_id
    parse_sample: str | None = None,
) -> pd.DataFrame:
    """
    Build a lookup table with one row per subject_id in this batch (Parse run),
    carrying the donor-level ASAP metadata (sample IDs, regions, sex) for joining
    to the AnnData object.

    SAMPLE.csv is one row per (sublibrary x sample), so it is first collapsed to one row
    per sample. A subject with more than one sample in the same batch (e.g. SN + ACG in
    PB-002) cannot be resolved from SAMPLE alone and needs sample_map, which records
    which Parse sample each ASAP sample was loaded in.
    """
    sample_cols = ["sample_id", "ASAP_sample_id", "subject_id", "batch",
                   "region_level_1", "region_level_2", "region_level_3"]

    # Collapse to one row per sample, asserting nothing varies across the sublibrary rows
    samples = sample_df[sample_cols].drop_duplicates()
    if not samples["sample_id"].is_unique:
        raise ValueError("Sample-level fields vary across pool_id rows in SAMPLE.csv")

    samples = samples[samples["batch"] == batch]
    if samples.empty:
        raise ValueError(f"No samples found in SAMPLE.csv for batch {batch}")

    if sample_map is not None:
        keep = sample_map.loc[sample_map["parse_sample"] == parse_sample, "sample_id"]
        ambiguous = samples["subject_id"].duplicated(keep=False)
        samples = samples[~ambiguous | samples["sample_id"].isin(keep)]

    dupes = samples.loc[samples["subject_id"].duplicated(), "subject_id"].unique()
    if len(dupes):
        raise ValueError(f"Subjects map to >1 sample in {batch}; need a sample map: {list(dupes)}")

    lookup = samples.drop(columns=["batch"]).merge(
        subject_df[["subject_id", "sex"]],
        on="subject_id",
        how="left",
        validate="one_to_one",
    )

    return lookup


def join_asap_metadata(
    adata: AnnData,
    sample_df: pd.DataFrame,
    subject_df: pd.DataFrame,
    batch: str,
    team: str,
    dataset: str,
    sample_map: pd.DataFrame | None = None,  # optional: parse_sample, subject_id, sample_id
    parse_sample: str | None = None,
) -> AnnData:
    """
    Join ASAP SAMPLE and SUBJECT metadata to the AnnData object, with field names and
    formatting matching sc-rnaseq-wf's prep_metadata:

    Run-level, assigned to every cell (including vireo doublet/unassigned):
    - batch, team, dataset, batch_id = f"{team}_{dataset}_{batch}"

    Donor-level, from the lookup (NA for vireo doublet/unassigned):
    - sample (ASAP_sample_id), sample_id, sex, brain_region_level_1/2/3

    sex and brain regions are lower-cased to match wf-common's generate_inputs.
    """
    lookup = build_sample_lookup(sample_df, subject_df, batch, sample_map, parse_sample)
    obs = adata.obs
    obs = obs.merge(lookup, on="subject_id", how="left", validate="many_to_one")
    obs.index = adata.obs.index  # merge drops the index; restore bc_wells

    # A missing sample is only expected for vireo's own non-singlet calls; anything 
    # else is a subject_id that does not match SAMPLE.csv for this batch
    vireo_non_singlets = {"doublet", "unassigned"}
    unmatched = obs["ASAP_sample_id"].isna() & ~obs["subject_id"].isin(vireo_non_singlets)
    if unmatched.any():
        missing = sorted(obs.loc[unmatched, "subject_id"].astype(str).unique())
        raise ValueError(
            f"{unmatched.sum()} cells have a subject_id not found in SAMPLE.csv for batch {batch}: {missing}"
        )
    
    obs = obs.rename(columns={
        "ASAP_sample_id": "sample",
        "region_level_1": "brain_region_level_1",
        "region_level_2": "brain_region_level_2",
        "region_level_3": "brain_region_level_3",
    })

    # Lower-case to match wf-common generate_inputs; tolerate missing values and all-NaN column
    for col in ["sex", "brain_region_level_1", "brain_region_level_2", "brain_region_level_3"]:
        obs[col] = obs[col].map(lambda x: x.lower() if isinstance(x, str) else x)

    # Run-level fields hold for every cell in the Parse run, so assign as scalars
    obs["batch"] = batch
    obs["team"] = team
    obs["dataset"] = dataset
    obs["batch_id"] = f"{team}_{dataset}_{batch}"

    n_non_singlet = obs["subject_id"].isin(vireo_non_singlets).sum()
    logging.info("Joined ASAP metadata; %d cells without a sample (vireo doublet/unassigned)", n_non_singlet)

    adata.obs = obs

    return adata


def main(args: argparse.Namespace):

    adata = create_adata_object(data_path=args.data_dir)
    adata = add_donor_assignments(adata, args.donor_assignments)
    adata = add_qc_metrics(adata)
    sample_df = pd.read_csv(args.sample_path)
    subject_df = pd.read_csv(args.subject_path)
    sample_map = pd.read_csv(args.sample_map) if args.sample_map else None
    parse_sample = str(adata.obs["parse_sample"].iloc[0])

    adata = join_asap_metadata(
        adata,
        sample_df=sample_df,
        subject_df=subject_df,
        batch=args.batch,
        team=args.team,
        dataset=args.dataset,
        sample_map=sample_map,
        parse_sample=parse_sample,
    )

    adata.write_h5ad(filename=args.output, compression="lzf")
    logging.info("Saved AnnData object to %s", args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build an AnnData object for one Parse sample from split-pipe "
                    "DGE output: attach donor assignments, compute QC metrics and "
                    "scrublet scores (as in sc-rnaseq-wf prep_metadata), and join "
                    "ASAP SAMPLE/SUBJECT metadata."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help="split-pipe DGE directory (count_matrix.mtx, all_genes.csv, cell_metadata.csv)"
    )
    parser.add_argument(
        "--donor-assignments",
        type=Path,
        required=True,
        help="Cell metadata CSV with subject_id from assign_donors.py"
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output h5ad"
    )
    parser.add_argument(
        "--sample-path",
        type=Path,
        required=True,
        help="ASAP SAMPLE.csv"
    )
    parser.add_argument(
        "--subject-path",
        type=Path,
        required=True,
        help="ASAP SUBJECT.csv"
    )
    parser.add_argument(
        "--batch",
        required=True,
        help="SAMPLE.csv batch value for this Parse run (e.g. Batch1)"
    )
    parser.add_argument(
        "--team",
        required=True,
        help="ASAP team ID; used in batch_id"
    )
    parser.add_argument(
        "--dataset",
        required=True,
        help="Dataset identifier (e.g. team-voet-pmdbs-sn-rnaseq-parsebio)"
    )
    parser.add_argument(
        "--sample-map",
        type=Path,
        help="CSV of parse_sample, subject_id, sample_id; required when a subject has "
             "more than one sample in the batch (e.g. PB-002)"
    )

    args = parser.parse_args()
    main(args)