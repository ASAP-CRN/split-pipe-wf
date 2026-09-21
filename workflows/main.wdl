version 1.0

# split-pipe scRNA-seq workflow entrypoint

import "structs.wdl"
import "split_pipe/split_pipe.wdl" as SplitPipe

workflow split_pipe_entrypoint {
	input {
		Array[Sublibrary] sublibraries
        String organism
        
        # Make reference genome if not provided
        File? genome_ref_tar_gz
        File genome_fasta
        File genes_gtf

        # Run split-pipe on the sublibraries and combine results
        String chemistry
        String kit
        File sample_list_tsv
        File? parameter_file
        
        String container_registry
        String zones = "us-central1-c us-central1-f"
	}
    
    call SplitPipe.split_pipe_wf {
        input:
            sublibraries = sublibraries,
            organism = organism,
            genome_ref_tar_gz = genome_ref_tar_gz,
            genome_fasta = genome_fasta,
            genes_gtf = genes_gtf,
            chemistry = chemistry,
            kit = kit,
            sample_list_tsv = sample_list_tsv,
            parameter_file = parameter_file,
            container_registry = container_registry,
            zones = zones
    }
    
    output {
        File combined_sublibs_tar_gz = split_pipe_wf.combined_sublibs_tar_gz
    }
    
    meta {
        description: "Main workflow to run split-pipe in the ASAP CRN Cloud environment."
    }

    parameter_meta {
        sublibraries: {help: "Array of sublibrary structs containing sublibrary IDs and FASTQ file paths."}
        organism: {help: "Organism; used to select which reference genome is used. Either 'human' (hg38) or 'mouse' (mm10)."}
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
