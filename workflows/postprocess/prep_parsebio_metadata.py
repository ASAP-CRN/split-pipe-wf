#!/usr/bin/env python3
# TODO Direct reference: https://www.sc-best-practices.org/preprocessing-visualization/quality-control/
# TODO: Doc and point to sc-rnaseq-wf

import argparse
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
    
    adata.var = gene_data.set_index("gene_name", drop=False)
    adata.var_names_make_unique()
    adata.var.index.name = "gene_symbol_unique"
    
    adata.obs = cell_meta_df.set_index("bc_wells", drop=False)

    logging.info("Created AnnData object with shape %s", adata.shape)
    return adata


# TODO: revist this logic w.r.t. join_donor_ids_from_vireo (also name sucks)
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

    vireo_meta_df = pd.read_csv(vireo_joined_metadata_path, index_col="bc_wells")
    if not adata.obs_names.equals(vireo_meta_df.index):
        raise ValueError("Cell order in the joined metadata does not match count_matrix.mtx rows")
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
    Build a lookup table with one row per subject_id in this batch (kit), carrying the
    donor-level ASAP metadata (sample IDs, regions, sex) for joining to the AnnData object.

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

    Kit-level, assigned to every cell (including vireo doublet/unassigned):
    - batch, team, dataset, batch_id = f"{team}_{dataset}_{batch}"

    Donor-level, from the lookup (NA for vireo doublet/unassigned):
    - sample (ASAP_sample_id), sample_id, sex, brain_region_level_1/2/3

    sex and brain regions are lower-cased to match wf-common's generate_inputs.
    """
    # split-pipe's "sample" is the Parse sample (e.g. pool_1); free the name for ASAP's sample
    obs = adata.obs.rename(columns={"sample": "parse_sample"})

    lookup = build_sample_lookup(sample_df, subject_df, batch, sample_map, parse_sample)
    obs = obs.merge(lookup, on="subject_id", how="left", validate="many_to_one")
    obs.index = adata.obs.index  # merge drops the index; restore bc_wells

    obs = obs.rename(columns={
        "ASAP_sample_id": "sample",
        "region_level_1": "brain_region_level_1",
        "region_level_2": "brain_region_level_2",
        "region_level_3": "brain_region_level_3",
    })

    # Lower-case to match wf-common generate_inputs; tolerate missing values and all-NaN column
    for col in ["sex", "brain_region_level_1", "brain_region_level_2", "brain_region_level_3"]:
        obs[col] = obs[col].map(lambda x: x.lower() if isinstance(x, str) else x)

    # Kit-level fields hold for every cell in the pool, so assign as scalars
    obs["batch"] = batch
    obs["team"] = team
    obs["dataset"] = dataset
    obs["batch_id"] = f"{team}_{dataset}_{batch}"

    n_unmatched = obs["sample"].isna().sum()
    logging.info("Joined ASAP metadata; %d cells without a sample (vireo doublet/unassigned)", n_unmatched)

    adata.obs = obs

    return adata


def main(args: argparse.Namespace):
    
    adata = create_adata_object(data_path=args.data_dir)
    adata = add_vireo_joined_metadata(adata, args.vireo_metadata_path)
    adata = add_qc_metrics(adata)
    sample_df = pd.read_csv(args.sample_path)
    subject_df = pd.read_csv(args.subject_path)
    
    adata = join_asap_metadata(
        adata, 
        sample_df=sample_df,
        subject_df=subject_df, 
        batch=args.batch, 
        team=args.team, 
        dataset=args.dataset,
        sample_map=None,  # TODO: Needs to be tested on PB-002
        parse_sample=None
    )
    
    adata.write_h5ad(filename=args.output, compression="gzip")
    logging.info("Saved AnnData object to %s", args.output)
    

if __name__ == "__main__":
    # TODO: doc 
    parser = argparse.ArgumentParser(description="TODO")
    
    parser.add_argument(
        "--data_dir", 
        type=Path,
        required=True, 
        help="Path to the directory containing the split-pipe output"
    )    
    parser.add_argument(
        "--vireo_metadata_path", 
        type=Path, 
        required=True,
        help="Path to the Vireo joined metadata CSV file"
    )    
    parser.add_argument(
        "--output",
        required=True,
        type=Path, 
        help="Path to the output object"
    )
    parser.add_argument(
        "--sample_path", 
        type=Path, 
        required=True,
        help="Path to SAMPLE.csv"
    )
    parser.add_argument(
        "--subject_path", 
        type=Path, 
        required=True,
        help="Path to SUBJECT.csv"
    )
    parser.add_argument(
        "--batch",
        required=True,
        help="Batch identifier for the samples"
    )
    parser.add_argument(
        "--team",
        required=True,
        help="Team identifier for the samples"
    )
    parser.add_argument(
        "--dataset",
        required=True,
        help="Dataset identifier for the samples"
    )
    
    args = parser.parse_args()
    main(args)