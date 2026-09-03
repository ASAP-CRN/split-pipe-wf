version 1.0

# Call split-pipe to generate a matrix of counts from paired-end FASTQ files

import "../structs.wdl"


workflow split_pipe_wf {
    input {
        Array[Sample] samples
        String organism
        File genome_fasta
        File genes_gtf

        String mode
        String chemistry
        String kit
        File sample_table
        File? parameter_file
		
        String container_registry
		String zones
    }

    Array[String] sublibrary_ids = samples.sublibrary_id

    # Build reference genome index
    call make_reference {
        input:
            organism = organism,
            genome_fasta = genome_fasta,
            genes_gtf = genes_gtf,
            container_registry = container_registry,
            zones = zones
    }

    # Run split-pipe on each sublibrary
    scatter (sample in samples) {
        call split_pipe {
            input:
                genome_ref_tar_gz = make_reference.genome_ref_tar_gz,
                mode = mode,
                sublibrary_id = sample.sublibrary_id,
                chemistry = chemistry,
                kit = kit,
                fastq_R1s = sample.fastq_R1s,
                fastq_R2s = sample.fastq_R2s,
                sample_table = sample_table,
                parameter_file = parameter_file
        }

    # Combine the output sublibaries into one object
    call combine_sublibraries {
        input:
            sublibrary_tars = split_pipe.sublibrary_tar_gz  # Gathered Array[File] from the split-pipe scatter
            sublibrary_ids = sublibrary_ids
    }

    output {
        # TODO: Output uncompressed directory with all outputs declared?
        File combined_sublibrary_tar_gz = combine_sublibraries.combined_sublibrary_tar_gz
    }

    meta {
        description = "This workflow create the indexed reference genome, runs split-pipe on each sublibrary and combines the results into a single object."
    }

    parameter_meta {
        samples: {help: "Array of sample structs containing sublibrary IDs and FASTQ file paths."}
        organism: {help: "Organism; used to select which reference genome is used. Either 'human' (hg38) or 'mouse' (mm10)."}
        genome_fasta: {help: "Path to the reference genome FASTA file."}
        genes_gtf: {help: "Path to the reference gene annotation GTF file."}
        mode: {help: "Mode to run split-pipe in."}
        chemistry: {help: "Chemistry used in the experiment."}
        kit: {help: "Kit used in the experiment."}
        sample_table: {help: "Path to the sample table file."}
        parameter_file: {help: "Optional path to the parameter file."}
        container_registry: {help: "Container registry where workflow Docker images are hosted."}
        zones: {help: "Compute zones to run the workflow in."}
    }
}


task make_reference {
    input {
        String organism
        File genome_fasta
        File genes_gtf
		String container_registry
		String zones
    }

    String genome_name = if organism == "human" then "hg38" else "mm10"

    Int threads = 16
    Int mem_gb = 48
    Int disk_size = 50

    command <<<
        set -euo pipefail

        mkdir genome_ref
        
        split-pipe \
            --mode mkref \
            --genome_name ~{genome_name} \
            --fasta ~{genome_fasta} \
            --genes ~{genes_gtf} \
            --output_dir genome_ref

        tar -czf genome_ref.tar.gz genome_ref
    >>>

    output {
        File genome_ref_tar_gz = "genome_ref.tar.gz"

    }

    runtime {
        docker: "~{container_registry}/split-pipe:1.8.2"
        cpu: threads
        cpuPlatform: "Intel Cascade Lake"  # TODO: necessary?
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} HDD"
        bootDiskSizeGb: 40
        zones: zones
    }

    meta {
        description = "Make reference genome index from FASTA and GTF files."
    }

    parameter_meta {
        organism: {help: "Organism; used to select which reference genome is used. Either 'human' (hg38) or 'mouse' (mm10)."}
        genome_fasta: {help: "Path to the reference genome FASTA file."}
        genes_gtf: {help: "Path to the reference gene annotation GTF file."}
    }

}


task split_pipe {
    input {
        File genome_ref_tar_gz = "genome_ref.tar.gz"   # TODO: outside of this block?
        String mode = "all"
        String sublibrary_id
        String chemistry
        String kit
        Array[File] fastq_R1s
        Array[File] fastq_R2s
        File sample_table
        File? parameter_file
        String container_registry
        String zones
    }

    Int threads = 16
    Int mem_gb = 48     # TODO: appropriate?
    Int disk_size = 50  # TODO: Much larger 

    command <<<
        set -euo pipefail

        # Extract the genome reference dir and make the output sublibrary dir
        tar -xzf ~{genome_ref_tar_gz}
        mkdir ~{sublibrary_id}

        split-pipe \
            --mode ~{mode} \
            --chemistry ~{chemistry} \
            --kit ~{kit} \
            --fq1 ~{fastq_R1s} \
            --fq2 ~{fastq_R2s} \
            --genome_dir genome_ref \
            --output_dir ~{sublibrary_id} \
            --samp_list ~{sample_table} \
            --parfile ~{parameter_file} \
            --dryrun      # To be removed after testing

        tar -czf ~{sublibrary_id}.tar.gz ~{sublibrary_id}        
    >>>

    output {
        File sublibrary_tar_gz = "~{sublibrary_id}.tar.gz"
    }

    runtime {
        docker: "~{container_registry}/split-pipe:1.8.2"
        cpu: threads
        cpuPlatform: "Intel Cascade Lake"  # TODO: necessary?
        memory: "{~mem_gb} GB"
        disks: "local-disk ~{disk_size} HDD"
        bootDiskSizeGb: 40
        zones: zones
    }

    meta {
        description = "Run split-pipe on the given sublibrary."  # TODO: more detail
    }

    # TODO: add options
    parameter_meta {
        sublibrary_id: {help: "Identifier for the sublibrary being processed."}
        chemistry: {help: "Chemistry used for the sublibrary."}
        kit: {help: "Kit used for the sublibrary."}
        fastqs_R1s: {help: "Array of FASTQ files for read 1."}
        fastqs_R2s: {help: "Array of FASTQ files for read 2."}
        sample_table: {help: "Path to the sample table file."}
        parameter_file: {help: "Optional path to the parameter file."}
    }
}


task combine_sublibraries {
    input {
        Array[File] sublibrary_tars
        Array[String] sublibrary_ids
        String container_registry
        String zones
    }
    
    # File of sublib dir paths to be fed as input to split-pipe combine
    File sublib_list = write_lines(sublibrary_ids)
    
    Int threads = 16
    Int mem_gb = 48     # TODO: appropriate?
    Int disk_size = 50  # TODO: Much larger

    command <<<
        set -euo pipefail

        # Extract each sublibrary .tar.gz into sublibs/
        mkdir -p sublibs
        ids=(~{sep=" " sublibrary_ids})  #  TODO: ask vs. write_lines()
        tars=(~{sep=" " sublibrary_tars})
        for i in "${!tars_array[@]}"; do
            mkdir -p "sublibs/${ids[$i]}"
            tar -xzf "${tars[$i]}" -C "sublibs/${ids[$i]}"
        done

        split-pipe \
            --mode comb \
            --sublib_list ~{sublib_list} \
            --sublib_pref "$PWD/sublibs/" \  # Prefix for the parent sublib dir
            --output_dir combined_sublibs

        # TODO: Full directory with all outputs declared?
        tar -czf combined_sublibs.tar.gz combined_sublibs
    >>>

    output {
        File combined_sublibs_tar_gz = "combined_sublibs.tar.gz"
    }

    runtime {
        docker: "~{container_registry}/split-pipe:1.8.2"
        cpu: threads
        cpuPlatform: "Intel Cascade Lake"  # TODO: necessary?
        memory: "{~mem_gb} GB"
        disks: "local-disk ~{disk_size} HDD"
        bootDiskSizeGb: 40
        zones: zones
    }

    meta {
        description = "Combine multiple sublibraries into a one combined object."
    }

    parameter_meta {
        sublibrary_tars: {help: "Array of sublibrary tar.gz files to be combined."}
        sublibrary_ids: {help: "Array of sublibrary identifiers corresponding to the tar.gz files."}
    }
