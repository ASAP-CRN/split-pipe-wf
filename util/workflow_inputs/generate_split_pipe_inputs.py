#!/usr/bin/env python3

"""
Generate the input JSON files required for running the split-pipe workflow, based
on the provided SAMPLE.csv, DATA.csv, and parse samples TSV files. 

This is currently tailored to the single released Parse Bio dataset (team-voet-pmdbs-sn-rnaseq-parsebio):
- Two sequencing kits/runs were used (PB-001 and PB-002) - run identifiers must 
  be tracked for merging into one final dataset.
- A Parse Bio run pools all samples/cells after barcoding - these are then split
  into sublibraries (identified by pool_id) for sequencing.
- The provided pool_id are of the form PB001a, where PB-001 is the run and PB-001a 
  is the sublibrary.
- This dataset also pools donors, requiring vireo de-multiplexing.
- These donor pools correspond to a "Parse sample" - what was initially aliquoted
  into the plate wells before barcoding/pooling.
- The exception is ASA_063B, which was not pooled and thus is its own "Parse sample."
- The sublibaries in the input *must* be specified in alphabetical order to be
  consistent with the provided vireo tables.

- Example:
    - Run: PB-001 [contains 4 Parse samples encompassing 16 donors,  split into 8 sublibraries]
    - Parse sample: pool_1 [contains donors ASA_004, ASA_005, ASA_006, ASA_041, ASA_042]
    - Sublibrary: PB-001a [FASTQs contain a proportion of cells from ALL donors in PB-001]
    
- This requires multiple tables to track these relationships:
    - Parse sample list (PB-001_sample_list.tsv): this tracks which Parse samples are
      present in each run's plate - necessary for split-pipe to resolve barcoding.
    - Vireo files (PB-001-pool_1.tsv): This associates each barcoded cell from a 
      Parse sample to its constituent donors, enabling donor de-multiplexing.
    - SAMPLE, SUBJECT, DATA: ASAP metadata tables that link samples to their FASTQs.
    - Sample map (PB-002_sample_map.csv): When a subject contributes more than one
      sample within a run (as in PB-002), this is needed to resolve the mapping
      from vireo since vireo only provides the donor/subject_id (i.e, ASA_004 not
      ASA_004A).
    - split_pipe_runs.tsv: used by this script to track Parse sample lists per run.
    - split_pipe_parse_samples.tsv: used by this script to associate vireo files
      (or the single donor ASA_063) to each Parse sample within a run.
"""

from datetime import datetime
from argparse import ArgumentParser
import pandas as pd
import re
import json


today = datetime.today()
date = today.strftime("%Y-%m-%d")

WORKFLOW_NAME = "split_pipe_analysis"

# Getting split-pipe --kit and --chemistry from SAMPLE.csv protocol_id
PROTOCOL_KIT_CHEMISTRY = {
    "Evercode_WT_v1": ("WT", "v1"),
    "Evercode_WT_MEGA_v1": ("WT_mega", "v1"),
}

# pool_id is the sublibrary, e.g. PB001a -> run PB-001, sublibrary PB-001a
POOL_ID_PATTERN = re.compile(r"^([A-Za-z]+)(\d+)([a-z]+)$")


def parse_pool_id(pool_id):
    match = POOL_ID_PATTERN.match(pool_id)
    if match is None:
        raise SystemExit(
            f"Unexpected pool_id format: '{pool_id}'. Expected <prefix><run number><sublibrary letter(s)>, e.g. PB001a."
        )
    prefix, run_number, sublibrary_letter = match.groups()
    run_id = f"{prefix}-{run_number}"
    sublibrary_id = f"{run_id}{sublibrary_letter}"
    return run_id, sublibrary_id


def pair_fastqs(pool_id, pool_df):
    """
    Return R1 and R2 FASTQ URIs for one sublibrary, sorted and paired by filename.
    DATA.csv repeats each file once per sample, so files are deduplicated first.
    """
    files = pool_df[["file_description", "gcp_uri"]].drop_duplicates()
    fastq_R1s = sorted(files.loc[files["file_description"].str.contains("R1"), "gcp_uri"])
    fastq_R2s = sorted(files.loc[files["file_description"].str.contains("R2"), "gcp_uri"])

    if len(fastq_R1s) == 0 or len(fastq_R1s) != len(fastq_R2s):
        raise SystemExit(
            f"Sublibrary {pool_id}: expected matching, non-zero R1 and R2 counts; found {len(fastq_R1s)} R1 and {len(fastq_R2s)} R2."
        )
    for fastq_R1, fastq_R2 in zip(fastq_R1s, fastq_R2s):
        # Only the last _R1_ is the read designator (sample names can contain R1 elsewhere)
        expected_R2 = "_R2_".join(fastq_R1.rsplit("_R1_", 1))
        if fastq_R2 != expected_R2:
            raise SystemExit(f"Sublibrary {pool_id}: R1 {fastq_R1} does not pair with R2 {fastq_R2}.")

    return fastq_R1s, fastq_R2s


def get_kit_chemistry(sample_df):
    """Return {batch: (kit, chemistry)} from SAMPLE.csv protocol_id; one protocol per batch."""
    kit_chemistry = dict()
    for batch, batch_df in sample_df.groupby("batch"):
        protocol_ids = sorted(batch_df["protocol_id"].dropna().unique())
        if len(protocol_ids) != 1:
            raise SystemExit(f"Batch {batch}: expected one protocol_id in SAMPLE.csv, found {protocol_ids}.")
        protocol_id = protocol_ids[0]
        if protocol_id not in PROTOCOL_KIT_CHEMISTRY:
            raise SystemExit(
                f"Batch {batch}: unknown protocol_id '{protocol_id}'; add it to PROTOCOL_KIT_CHEMISTRY. "
                f"Known: {sorted(PROTOCOL_KIT_CHEMISTRY)}"
            )
        kit_chemistry[batch] = PROTOCOL_KIT_CHEMISTRY[protocol_id]
    return kit_chemistry


def get_parse_samples(parse_samples_df, run_id):
    """Parse samples for one run, in TSV order; each has exactly one of vireo_donor_ids or subject_id."""
    run_df = parse_samples_df[parse_samples_df["run_id"] == run_id]
    if run_df.empty:
        raise SystemExit(f"Run {run_id}: no rows in the parse samples TSV.")
    if run_df["parse_sample_name"].duplicated().any():
        raise SystemExit(f"Run {run_id}: duplicate parse_sample_name in the parse samples TSV.")

    parse_samples = []
    for row in run_df.itertuples():
        has_vireo = not pd.isna(row.vireo_donor_ids)
        has_subject = not pd.isna(row.subject_id)
        if has_vireo == has_subject:
            raise SystemExit(
                f"Run {run_id}, Parse sample {row.parse_sample_name}: set exactly one of vireo_donor_ids (pooled) or subject_id (un-pooled)."
            )
        parse_sample = {"parse_sample_name": str(row.parse_sample_name)}
        if has_vireo:
            parse_sample["vireo_donor_ids"] = row.vireo_donor_ids
        else:
            parse_sample["subject_id"] = row.subject_id
        parse_samples.append(parse_sample)
    return parse_samples


def main(args):
    data_df = pd.read_csv(args.data_csv)
    # ASAP dataset ID (e.g. DS_PMDBS_0038) is the workflow dataset_id, as in sc-rnaseq-wf;
    # --dataset-id (the bucket-style name) is only used for bucket paths and the output filename
    asap_dataset_ids = sorted(data_df["ASAP_dataset_id"].unique())
    if len(asap_dataset_ids) != 1:
        raise SystemExit(f"Expected one ASAP_dataset_id in DATA.csv, found {asap_dataset_ids}.")
    asap_dataset_id = asap_dataset_ids[0]
    sample_df = pd.read_csv(args.sample_csv)
    runs_df = pd.read_csv(args.runs_tsv, delimiter="\t", dtype=str)
    parse_samples_df = pd.read_csv(args.parse_samples_tsv, delimiter="\t", dtype=str)

    kit_chemistry = get_kit_chemistry(sample_df)

    data_df = data_df[data_df["file_type"] == "fastq"].copy()
    data_df[["run_id", "sublibrary_id"]] = data_df["pool_id"].apply(lambda pool_id: pd.Series(parse_pool_id(pool_id)))

    data_run_ids = sorted(data_df["run_id"].unique())
    config_run_ids = sorted(runs_df["run_id"].unique())
    if data_run_ids != config_run_ids:
        raise SystemExit(f"Runs in DATA.csv {data_run_ids} do not match runs in the runs TSV {config_run_ids}.")

    parse_runs = []
    for run_id in data_run_ids:
        run_df = data_df[data_df["run_id"] == run_id]

        batches = sorted(run_df["batch"].unique())
        if len(batches) != 1:
            raise SystemExit(f"Run {run_id}: expected one batch in DATA.csv, found {batches}.")
        batch = batches[0]
        if batch not in kit_chemistry:
            raise SystemExit(f"Run {run_id}: batch {batch} not found in SAMPLE.csv.")
        kit, chemistry = kit_chemistry[batch]

        # Sorted alphabetically: comb's __s<n> cell suffixes are positional and the vireo
        # tables were computed on sublibraries in alphabetical order
        sublibraries = []
        for pool_id, pool_df in sorted(run_df.groupby("pool_id"), key=lambda item: item[0]):
            fastq_R1s, fastq_R2s = pair_fastqs(pool_id, pool_df)
            sublibraries.append({
                "sublibrary_id": pool_df["sublibrary_id"].iloc[0],
                "fastq_R1s": fastq_R1s,
                "fastq_R2s": fastq_R2s,
            })

        run_config = runs_df[runs_df["run_id"] == run_id].iloc[0]
        parse_run = {
            "run_id": run_id,
            "kit": kit,
            "chemistry": chemistry,
            "batch": batch,
            "sample_list_tsv": run_config["sample_list_tsv"],
        }
        if not pd.isna(run_config.get("sample_map")):
            parse_run["sample_map"] = run_config["sample_map"]
        parse_run["sublibraries"] = sublibraries
        parse_run["parse_samples"] = get_parse_samples(parse_samples_df, run_id)
        parse_runs.append(parse_run)

    with open(args.inputs_template, "r") as f:
        inputs_json = json.load(f)

    staging_env = "dev" if args.embargoed else "uat"
    staging_bucket_types = ["dev"] if args.embargoed else ["dev", "uat"]

    inputs_json[f"{WORKFLOW_NAME}.dataset_id"] = asap_dataset_id
    inputs_json[f"{WORKFLOW_NAME}.team_id"] = args.team_id
    inputs_json[f"{WORKFLOW_NAME}.parse_runs"] = parse_runs
    inputs_json[f"{WORKFLOW_NAME}.raw_data_bucket"] = f"gs://asap-raw-{args.dataset_id}"
    inputs_json[f"{WORKFLOW_NAME}.staging_data_buckets"] = [
        f"gs://asap-{staging_bucket_type}-{args.dataset_id}"
        for staging_bucket_type in staging_bucket_types
    ]

    # Catch template entries that are still type descriptions rather than values
    unfilled = [
        key for key, value in inputs_json.items()
        if isinstance(value, str) and re.match(r"^(String|File|Int|Float|Boolean|Array)\b", value)
    ]
    if unfilled:
        print(f"WARNING: template values still look like type descriptions: {unfilled}")

    inputs_output_filename = f"inputs.{staging_env}.{args.dataset_id}.{date}.json"
    with open(inputs_output_filename, "w") as f:
        json.dump(inputs_json, f, indent=4)
    print(f"Wrote input JSON file: {inputs_output_filename}")
    print(f"  dataset_id: {asap_dataset_id}")
    for parse_run in parse_runs:
        print(
            f"  {parse_run['run_id']}: {parse_run['batch']}, {parse_run['kit']}/{parse_run['chemistry']}, "
            f"{len(parse_run['sublibraries'])} sublibraries, {len(parse_run['parse_samples'])} Parse samples"
        )


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Given ASAP DATA.csv and SAMPLE.csv plus per-run Parse sample information, generate a split-pipe-wf inputs JSON"
    )

    parser.add_argument(
        "-d",
        "--data-csv",
        type=str,
        required=True,
        help="ASAP DATA.csv; FASTQ gs URIs per sublibrary (pool_id), batch per sublibrary.",
    )
    parser.add_argument(
        "-s",
        "--sample-csv",
        type=str,
        required=True,
        help="ASAP SAMPLE.csv; protocol_id per batch gives the split-pipe kit and chemistry.",
    )
    parser.add_argument(
        "-r",
        "--runs-tsv",
        type=str,
        required=True,
        help="One row per Parse run; columns run_id, sample_list_tsv, sample_map (optional, may be empty).",
    )
    parser.add_argument(
        "-p",
        "--parse-samples-tsv",
        type=str,
        required=True,
        help="One row per Parse sample; columns run_id, parse_sample_name, vireo_donor_ids, subject_id (set exactly one of the last two).",
    )
    parser.add_argument(
        "-i",
        "--inputs-template",
        type=str,
        required=True,
        help="Inputs JSON with run-independent values filled in (reference, MMC, filtering parameters, container_registry); dataset, run and bucket keys are overwritten.",
    )
    parser.add_argument(
        "-t",
        "--dataset-id",
        type=str,
        required=True,
        help="Dataset name used in bucket names (gs://asap-raw-{dataset_id}, gs://asap-{dev,uat}-{dataset_id}), e.g. team-voet-pmdbs-sn-rnaseq-parsebio.",
    )
    parser.add_argument(
        "-a",
        "--team-id",
        type=str,
        required=True,
        help="ASAP team identifier; used in batch_id.",
    )
    parser.add_argument(
        "-e",
        "--embargoed",
        action="store_true",
        required=False,
        help="Dataset is embargoed: stage to the dev bucket only. Otherwise stage to dev and uat [false].",
    )

    args = parser.parse_args()

    main(args)
