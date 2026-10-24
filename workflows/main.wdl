version 1.0

# split-pipe scRNA-seq workflow entrypoint

import "structs.wdl"
import "split_pipe/split_pipe.wdl" as SplitPipe
import "split_pipe/postprocess.wdl" as PostProcess

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

        Boolean merge_runs = true
        
        String container_registry
        String zones = "us-central1-c us-central1-f"
    }
    
    # Build once outside of the per-run scatter
    if (!defined(genome_ref_tar_gz)) {
        call SplitPipe.make_genome_ref {
            input:
                genome_name = genome_name,
                genome_fasta = genome_fasta,
                genes_gtf = genes_gtf,
                container_registry = container_registry,
                zones = zones
        }
    }

    File genome_ref = select_first([genome_ref_tar_gz, make_genome_ref.genome_ref_tar_gz])
    
    scatter run in parse_runs {
        call SplitPipe.split_pipe {
            input:
                dataset_id = dataset_id,
                sublibraries = run.sublibraries,
                genome_name = genome_name,
                genome_ref_tar_gz = genome_ref_tar_gz,
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

    if (merge_runs) {
        call merge_runs {
            input:
                cleaned_run_adata_objects = postprocess.cleaned_adata,
                dataset_id = dataset_id,
                team_id = team_id,
                container_registry = container_registry,
                zones = zones
        }
    }

     output {
        File final_merged_adata = postprocess.merged_adata
    }

    meta {
        description: "Main workflow to run split-pipe in the ASAP CRN Cloud environment."
    }

    parameter_meta {
        dataset_id: {help: "The unique identifier for the dataset."}
        sublibraries: {help: "Array of sublibrary structs containing sublibrary IDs and FASTQ file paths."}
        genome_name: {help: "The name of which reference genome is used. Either 'hg38' (human) or 'mm10' (mouse)."}
        genome_fasta: {help: "Path to the reference genome FASTA file."}
        genes_gtf: {help: "Path to the reference gene annotation GTF file."}
        chemistry: {help: "Chemistry used in the experiment."}
        kit: {help: "Kit used in the experiment."}
        sample_list_tsv: {help: "Path to the sample loading table file."}
        parameter_file: {help: "Optional path to the parameter file."}
        container_registry: {help: "Container registry where workflow Docker images are hosted."}
        zones: {help: "Compute zones to run the workflow in."}
    }

}


task merge_runs {
    input {
        Array[File] cleaned_run_adata_objects
        String dataset_id
        String team_id
        String container_registry
        String zones
    }

    command {
        echo "Merging runs with output prefix: ${output_prefix}"
    }
    
    output {
        File merged_adata_object = "${dataset_id}_final_adata.h5ad"
    }
}