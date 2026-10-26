version 1.0

# split-pipe scRNA-seq workflow entrypoint

import "structs.wdl"
import "split_pipe/split_pipe.wdl" as SplitPipe
import "postprocess/postprocess.wdl" as PostProcess

workflow split_pipe_entrypoint {
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
        
        String container_registry
        String zones = "us-central1-c us-central1-f"
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
                sublibraries = run.sublibraries,
                genome_name = genome_name,
                genome_ref_tar_gz = genome_ref,
                genome_fasta = genome_fasta,
                genes_gtf = genes_gtf,
                chemistry = run.chemistry,
                kit = run.kit,
                sample_list_tsv = run.sample_list_tsv,
                parameter_file = parameter_file,
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
                container_registry = container_registry,
                zones = zones
        }
    }

    call merge_runs {
        input:
            cleaned_run_adata_objects = postprocess.mmc_adata_object,
            dataset_id = dataset_id,
            container_registry = container_registry,
            zones = zones
    }

    output {
        File final_merged_adata = merge_runs.merged_adata_object
    }

    meta {
        description: "Main workflow to run split-pipe in the ASAP CRN Cloud environment."
    }

    parameter_meta {
        dataset_id: {help: "The unique identifier for the dataset."}
        team_id: {help: "ASAP team identifier"}
        parse_runs: {help: "Array of ParseRun structs containing run-level information for each sequencing run."}
        genome_name: {help: "The name of which reference genome is used. Either 'hg38' (human) or 'mm10' (mouse)."}
        genome_fasta: {help: "Path to the reference genome FASTA file."}
        genes_gtf: {help: "Path to the reference gene annotation GTF file."}
        genome_ref_tar_gz: {help: "Optional prebuilt split-pipe mkref index (tar.gz). If omitted, one is built from genome_fasta and genes_gtf."}
        parameter_file: {help: "Optional path to the parameter file."}
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
        container_registry: {help: "Container registry where workflow Docker images are hosted."}
        zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
    }
}


task merge_runs {
    input {
        Array[File] cleaned_run_adata_objects
        String dataset_id
        String container_registry
        String zones
    }

    Int mem_gb = ceil(size(cleaned_run_adata_objects, "GB") * 6 + 20)
    Int disk_size = ceil(size(cleaned_run_adata_objects, "GB") * 3 + 50)

    command <<<
        set -euo pipefail
        
        merge_runs.py \
            --adata-inputs ~{sep=' ' cleaned_run_adata_objects} \
            --adata-output ~{dataset_id}.merged.mmc.h5ad
    >>>
    
    output {
        File merged_adata_object = "~{dataset_id}.merged.mmc.h5ad"
    }

    runtime {
        docker: "~{container_registry}/parsebio_tools:0.3.0"
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
        container_registry: {help: "Container registry where workflow Docker images are hosted."}
        zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
    }
}