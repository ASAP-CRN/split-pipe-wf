version 1.0

# split-pipe scRNA-seq workflow entrypoint

import "structs.wdl"
import "../wf-common/wdl/tasks/get_workflow_metadata.wdl" as GetWorkflowMetadata
import "../wf-common/wdl/tasks/upload_final_outputs.wdl" as UploadFinalOutputs
import "split_pipe/split_pipe.wdl" as SplitPipe
import "postprocess/postprocess.wdl" as PostProcess

workflow split_pipe_analysis {
    input {
        String dataset_id
        String team_id
        Array[ParseRun] parse_runs

        # Make reference genome if not provided
        String genome_name
        File genome_fasta
        File genes_gtf
        File? genome_ref_tar_gz

        # Optional parameter file for configuring split-pipe run
        File? parameter_file
        
        # Tables to map cell IDs -> samples
        File sample_csv
        File subject_csv

        # Filtering parameters
        Int pct_counts_mt_max = 10
        Float doublet_score_max = 0.2
        Array[Int] total_counts_limits = [500, 150000]
        Array[Int] n_genes_by_counts_limits = [300, 15000]
        Float n_mads_lower = 3
        Float n_mads_upper = 5

        # Allen Institute's Map My Cells
        String mmc_taxonomy
        File allen_brain_mmc_precomputed_stats_h5
        File? allen_brain_mmc_marker_genes_json
        
        String raw_data_bucket
        Array[String] staging_data_buckets
        String container_registry
        String zones = "us-central1-c us-central1-f"
    }

    String workflow_name = "pmdbs_split_pipe"
    String workflow_execution_path = "workflow_execution"
    String workflow_version = "v1.0.0"
    String workflow_release = "https://github.com/ASAP-CRN/split-pipe-wf/releases/tag/split_pipe_analysis-~{workflow_version}"
    String crn_release_version = "v6.0.0"
    String raw_data_path_prefix = "~{raw_data_bucket}/~{workflow_execution_path}/~{workflow_name}"
    String staging_path_prefix = "~{workflow_name}/release/~{crn_release_version}"


    call GetWorkflowMetadata.get_workflow_metadata {
        input:
            zones = zones
    }

    # Build once outside of the per-run scatter so split_pipe skips its own make_reference
    if (!defined(genome_ref_tar_gz)) {
        call SplitPipe.make_reference {
            input:
                genome_name = genome_name,
                genome_fasta = genome_fasta,
                genes_gtf = genes_gtf,
                container_registry = container_registry,
                zones = zones
        }
    }

    File genome_ref = select_first([genome_ref_tar_gz, make_reference.genome_ref_tar_gz])
    
    scatter (run in parse_runs) {
        call SplitPipe.split_pipe {
            input:
                dataset_id = dataset_id,
                run_id = run.run_id,
                sublibraries = run.sublibraries,
                genome_name = genome_name,
                genome_ref_tar_gz = genome_ref,
                genome_fasta = genome_fasta,
                genes_gtf = genes_gtf,
                chemistry = run.chemistry,
                kit = run.kit,
                sample_list_tsv = run.sample_list_tsv,
                parameter_file = parameter_file,
                workflow_name = workflow_name,
                workflow_version = workflow_version,
                workflow_release = workflow_release,
                run_timestamp = get_workflow_metadata.timestamp,
                raw_data_path_prefix = raw_data_path_prefix,
                billing_project = get_workflow_metadata.billing_project,
                container_registry = container_registry,
                zones = zones
        }

        call PostProcess.postprocess {
            input:
                dataset_id = dataset_id,
                team_id = team_id,
                combined_sublibs_tar_gz = split_pipe.combined_sublibs_tar_gz,
                run_id = run.run_id,
                batch = run.batch,
                parse_samples = run.parse_samples,
                sample_csv = sample_csv,
                subject_csv = subject_csv,
                sample_map = run.sample_map,
                pct_counts_mt_max = pct_counts_mt_max,
                doublet_score_max = doublet_score_max,
                total_counts_limits = total_counts_limits,
                n_genes_by_counts_limits = n_genes_by_counts_limits,
                n_mads_lower = n_mads_lower,
                n_mads_upper = n_mads_upper,
                mmc_taxonomy = mmc_taxonomy,
                allen_brain_mmc_precomputed_stats_h5 = allen_brain_mmc_precomputed_stats_h5,
                allen_brain_mmc_marker_genes_json = allen_brain_mmc_marker_genes_json,
                workflow_name = workflow_name,
                workflow_version = workflow_version,
                workflow_release = workflow_release,
                run_timestamp = get_workflow_metadata.timestamp,
                raw_data_path_prefix = raw_data_path_prefix,
                billing_project = get_workflow_metadata.billing_project,
                container_registry = container_registry,
                zones = zones
        }
    }

    Array[Array[String]] workflow_info = [[get_workflow_metadata.timestamp, workflow_name, workflow_version, workflow_release]]
    # Put merged output in postprocess dir - each run has same path so just select the first
    String postprocess_raw_path = postprocess.postprocess_raw_data_path[0]

    call merge_runs {
        input:
            cleaned_run_adata_objects = postprocess.mmc_adata_object,
            dataset_id = dataset_id,
            raw_data_path = postprocess_raw_path,
            workflow_info = workflow_info,
            billing_project = get_workflow_metadata.billing_project,
            container_registry = container_registry,
            zones = zones
    }

    # Uploaded outputs into staging
    Array[String] split_pipe_final_output_paths = flatten([
        flatten(split_pipe.sublibrary_bams),
        flatten(split_pipe.sublibrary_star_logs),
        split_pipe.combined_sublibs_tar_gz,
        split_pipe.all_sample_summary_html,
        split_pipe.agg_sample_summary_csv
    ]) #!StringCoercion

    Array[String] postprocess_final_output_paths = flatten([
        postprocess.merged_adata,
        postprocess.initial_metadata,
        flatten(postprocess.qc_plots),
        postprocess.unassigned_removed_cells_csv,
        postprocess.filtered_out_cells_csv,
        postprocess.mmc_extended_results_json,
        postprocess.mmc_results_csv,
        postprocess.mmc_log_txt,
        postprocess.mmc_results_parquet,
        [merge_runs.merged_run_adata_object]
    ]) #!StringCoercion

    call UploadFinalOutputs.upload_final_outputs as upload_split_pipe_files {
        input:
            output_file_paths = split_pipe_final_output_paths,
            staging_data_buckets = staging_data_buckets,
            staging_data_path = "~{staging_path_prefix}/split_pipe",
            billing_project = get_workflow_metadata.billing_project,
            zones = zones
    }

    call UploadFinalOutputs.upload_final_outputs as upload_postprocess_files {
        input:
            output_file_paths = postprocess_final_output_paths,
            staging_data_buckets = staging_data_buckets,
            staging_data_path = "~{staging_path_prefix}/postprocess",
            billing_project = get_workflow_metadata.billing_project,
            zones = zones
    }

    output {
        # Reference - pass back as genome_ref_tar_gz to skip make_reference in later submissions
        File genome_ref_tar_gz_used = genome_ref
        # Split-pipe
        Array[Array[File]] sublibrary_tar_gz = split_pipe.sublibrary_tar_gz
        Array[Array[File]] sublibrary_bams = split_pipe.sublibrary_bams
        Array[Array[File]] sublibrary_star_logs = split_pipe.sublibrary_star_logs
        Array[File] combined_sublibs_tar_gz = split_pipe.combined_sublibs_tar_gz
        Array[File] all_sample_summary_html = split_pipe.all_sample_summary_html
        Array[File] agg_sample_summary_csv = split_pipe.agg_sample_summary_csv
        # Postprocess
        Array[Array[File]] prepped_adata = postprocess.prepped_adata
        Array[File] merged_adata = postprocess.merged_adata
        Array[File] initial_metadata = postprocess.initial_metadata
        Array[Array[File]] qc_plots = postprocess.qc_plots
        Array[File] filtered_adata_object = postprocess.filtered_adata_object
        Array[File] unassigned_removed_cells_csv = postprocess.unassigned_removed_cells_csv
        Array[File] filtered_out_cells_csv = postprocess.filtered_out_cells_csv
        Array[File] mmc_extended_results_json = postprocess.mmc_extended_results_json
        Array[File] mmc_results_csv = postprocess.mmc_results_csv
        Array[File] mmc_log_txt = postprocess.mmc_log_txt
        Array[File] mmc_adata_object = postprocess.mmc_adata_object
        Array[File] mmc_results_parquet = postprocess.mmc_results_parquet
        # Merged runs
        File final_merged_adata = merge_runs.merged_run_adata_object #!FileCoercion
        # Staging manifests
        Array[File] split_pipe_manifest_tsvs = upload_split_pipe_files.manifests #!FileCoercion
        Array[File] postprocess_manifest_tsvs = upload_postprocess_files.manifests #!FileCoercion
    }

    meta {
        description: "Main workflow to run split-pipe in the ASAP CRN Cloud environment."
    }

    parameter_meta {
        dataset_id: {help: "ASAP dataset ID (e.g. DS_PMDBS_0038); used in output names and the obs dataset and batch_id columns."}
        team_id: {help: "ASAP team identifier"}
        parse_runs: {help: "Array of ParseRun structs containing run-level information for each sequencing run."}
        genome_name: {help: "The name of which reference genome is used. Either 'hg38' (human) or 'mm10' (mouse)."}
        genome_fasta: {help: "Path to the reference genome FASTA file."}
        genes_gtf: {help: "Path to the reference gene annotation GTF file."}
        genome_ref_tar_gz: {help: "Optional prebuilt split-pipe mkref index (tar.gz). If omitted, one is built from genome_fasta and genes_gtf."}
        parameter_file: {help: "Optional split-pipe --parfile; 'key value' per line."}
        sample_csv: {help: "ASAP SAMPLE.csv."}
        subject_csv: {help: "ASAP SUBJECT.csv."}
        pct_counts_mt_max: {help: "Maximum percentage of mitochondrial gene counts allowed per cell. [10]"}
        doublet_score_max: {help: "Maximum doublet detection score threshold. [0.2]"}
        total_counts_limits: {help: "Absolute minimum and maximum total UMI counts per cell, applied on top of the MAD-based thresholds. [500, 150000]"}
        n_genes_by_counts_limits: {help: "Absolute minimum and maximum number of genes detected per cell, applied on top of the MAD-based thresholds. [300, 15000]"}
        n_mads_lower: {help: "Number of MADs below the per-sample median allowed for total counts and genes detected. [3]"}
        n_mads_upper: {help: "Number of MADs above the per-sample median allowed for total counts and genes detected. [5]"}
        mmc_taxonomy: {help: "Cell type taxonomy of the precomputed stats reference; appended to MMC output filenames. Must match allen_brain_mmc_precomputed_stats_h5."}
        allen_brain_mmc_precomputed_stats_h5: {help: "A precomputed statistics file from the Allen Brain Cell Atlas containing reference statistics (the average gene expression profile per cell type cluster and cell type taxonomy)."}
        allen_brain_mmc_marker_genes_json: {help: "Optional marker gene lookup (parent node -> marker genes), e.g. extracted from a previous run's extended_results.json. If provided, skips on-the-fly marker computation."}
        raw_data_bucket: {help: "Dataset raw data bucket (e.g. gs://asap-raw-team-voet-pmdbs-sn-rnaseq-parsebio); outputs are uploaded under `<raw_data_bucket>/workflow_execution/<workflow_name>/`."}
        staging_data_buckets: {help: "Set of staging buckets (dev and uat) for uploading workflow outputs."}
        container_registry: {help: "Container registry where workflow Docker images are hosted."}
        zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
    }
}


task merge_runs {
    input {
        Array[File] cleaned_run_adata_objects
        String dataset_id

        String raw_data_path
        Array[Array[String]] workflow_info
        String billing_project
        String container_registry
        String zones
    }

    String merged_run_adata_name = "~{dataset_id}.merged.mmc.h5ad"

    Int mem_gb = ceil(size(cleaned_run_adata_objects, "GB") * 6 + 20)
    Int disk_size = ceil(size(cleaned_run_adata_objects, "GB") * 3 + 50)

    command <<<
        set -euo pipefail
        
        /usr/bin/time -v \
        merge_runs.py \
            --adata-inputs ~{sep=' ' cleaned_run_adata_objects} \
            --adata-output ~{merged_run_adata_name}

        upload_outputs \
            -b ~{billing_project} \
            -d ~{raw_data_path} \
            -i ~{write_tsv(workflow_info)} \
            -o ~{merged_run_adata_name}
    >>>
    
    output {
        String merged_run_adata_object = "~{raw_data_path}/~{merged_run_adata_name}"
    }

    runtime {
        docker: "~{container_registry}/parsebio_tools:0.3.3"
        cpu: 4
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 3
        bootDiskSizeGb: 40
        zones: zones
    }

    meta {
        description: "Merges the cleaned run-level Parse Bio AnnData objects into a single AnnData object."
    }

    parameter_meta {
        cleaned_run_adata_objects: {help: "An array of cleaned run-level Parse Bio AnnData objects to merge."}
        dataset_id: {help: "The unique identifier for the dataset; used to name the output."}
        raw_data_path: {help: "Raw data bucket path for the merged run outputs."}
        workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
        billing_project: {help: "Billing project to charge GCP costs."}
        container_registry: {help: "Container registry where workflow Docker images are hosted."}
        zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
    }
}