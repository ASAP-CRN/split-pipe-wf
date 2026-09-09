version 1.0

# Call split-pipe to generate a matrix of counts from paired-end FASTQ files

import "../structs.wdl"


workflow split_pipe_wf {
    input {
        Array[Sublibrary] sublibraries
        String organism
        File genome_fasta
        File genes_gtf

        String chemistry
        String kit
        File sample_list_tsv
        File? parameter_file
		
        String container_registry
		String zones
    }

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
    scatter (sublibrary in sublibraries) {
        String sublibrary_id = sublibrary.sublibrary_id  # Declare as gathered for array input for combine_sublibraries
        
        call split_pipe {
            input:
                genome_ref_tar_gz = make_reference.genome_ref_tar_gz,
                sublibrary_id = sublibrary_id,
                chemistry = chemistry,
                kit = kit,
                fastq_R1s = sublibrary.fastq_R1s,
                fastq_R2s = sublibrary.fastq_R2s,
                sample_list_tsv = sample_list_tsv,
                parameter_file = parameter_file,
                container_registry = container_registry,
                zones = zones
        }
    }

    Array[String] sublibrary_ids = sublibrary_id

    # Combine the output sublibaries into one object
    call combine_sublibraries {
        input:
            sublibrary_tars = split_pipe.sublibrary_tar_gz,  # Gathered Array[File] from the split-pipe scatter
            sublibrary_ids = sublibrary_ids,
            container_registry = container_registry,
            zones = zones
    }

    output {
        # TODO: Output uncompressed directory with all outputs declared?
        File combined_sublibs_tar_gz = combine_sublibraries.combined_sublibs_tar_gz
    }

    meta {
        description: "This workflow create the indexed reference genome, runs split-pipe on each sublibrary and combines the results into a single object."
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
            --output_dir genome_ref \
            --nthreads ~{threads}

        tar -czf genome_ref.tar.gz genome_ref
    >>>

    output {
        File genome_ref_tar_gz = "genome_ref.tar.gz"

    }

    runtime {
        docker: "~{container_registry}/split_pipe:1.8.2"
        cpu: threads
        cpuPlatform: "Intel Cascade Lake"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} HDD"
        bootDiskSizeGb: 40
        zones: zones
    }

    meta {
        description: "Make reference genome index from FASTA and GTF files."
    }

    parameter_meta {
        organism: {help: "Organism; used to select which reference genome is used. Either 'human' (hg38) or 'mouse' (mm10)."}
        genome_fasta: {help: "Path to the reference genome FASTA file."}
        genes_gtf: {help: "Path to the reference gene annotation GTF file."}
    }

}


task split_pipe {
    input {
        File genome_ref_tar_gz
        String sublibrary_id
        String chemistry
        String kit
        Array[File] fastq_R1s
        Array[File] fastq_R2s
        File sample_list_tsv
        File? parameter_file
        String container_registry
        String zones
    }

    Int threads = 16
    Int mem_gb = 48
    Int disk_size = ceil((size(genome_ref_tar_gz, "GB") + size(flatten([fastq_R1s, fastq_R2s]), "GB")) * 4 + 50)

    command <<<
        set -euo pipefail

        # Extract the genome reference dir and make the output sublibrary dir
        tar -xzf ~{genome_ref_tar_gz}
        mkdir ~{sublibrary_id}

        split-pipe \
            --mode all \
            --chemistry ~{chemistry} \
            --kit ~{kit} \
            --fq1 ~{sep=" " fastq_R1s} \
            --fq2 ~{sep=" " fastq_R2s} \
            --genome_dir genome_ref \
            --output_dir ~{sublibrary_id} \
            --samp_list ~{sample_list_tsv} \
            --nthreads ~{threads} \
            --no_keep_going \
            ~{if defined(parameter_file) then "--parfile " + select_first([parameter_file]) else ""} \
            --dryrun      # TODO: remove after testing

        tar -czf ~{sublibrary_id}.tar.gz ~{sublibrary_id}        
    >>>

    output {
        File sublibrary_tar_gz = "~{sublibrary_id}.tar.gz"
    }

    runtime {
        docker: "~{container_registry}/split_pipe:1.8.2"
        cpu: threads
        cpuPlatform: "Intel Cascade Lake"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} HDD"
        bootDiskSizeGb: 40
        zones: zones
    }

    meta {
        description: "Run split-pipe on the given sublibrary."  # TODO: more detail
    }

    # TODO: add options
    parameter_meta {
        sublibrary_id: {help: "Identifier for the sublibrary being processed."}
        chemistry: {help: "Chemistry used for the sublibrary."}
        kit: {help: "Kit used for the sublibrary."}
        fastq_R1s: {help: "Array of FASTQ files for read 1."}
        fastq_R2s: {help: "Array of FASTQ files for read 2."}
        sample_list_tsv: {help: "Path to the sample table file."}
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
    
    File sublib_list = write_lines(sublibrary_ids)
    Int threads = 16
    Int mem_gb = 48  # TODO: Revisit, PB docs bin at 64/128/256 GB
    Int disk_size = ceil(size(sublibrary_tars, "GB") * 8 + 50)  # TODO: 8x for unpacking each tar + the final combined object - revisit 

    command <<<
        set -euo pipefail

        # Extract each sublibrary .tar.gz into sublibs/
        mkdir -p sublibs
        ids_array=(~{sep=" " sublibrary_ids})
        tars_array=(~{sep=" " sublibrary_tars})
        for i in "${!tars_array[@]}"; do
            mkdir -p "sublibs/${ids_array[$i]}"
            tar -xzf "${tars_array[$i]}" -C "sublibs/${ids_array[$i]}"
        done

        split-pipe \
            --mode comb \
            --sublib_list ~{sublib_list} \
            --sublib_pref "$PWD/sublibs/" \
            --output_dir combined_sublibs \
            --nthreads ~{threads}

        # TODO: Full directory with all outputs declared?
        tar -czf combined_sublibs.tar.gz combined_sublibs
    >>>

    output {
        File combined_sublibs_tar_gz = "combined_sublibs.tar.gz"
    }

    runtime {
        docker: "~{container_registry}/split_pipe:1.8.2"
        cpu: threads
        cpuPlatform: "Intel Cascade Lake"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} HDD"
        bootDiskSizeGb: 40
        zones: zones
    }

    meta {
        description: "Combine multiple sublibraries into a single combined object."
    }

    parameter_meta {
        sublibrary_tars: {help: "Array of sublibrary tar.gz files to be combined."}
        sublibrary_ids: {help: "Array of sublibrary identifiers corresponding to the tar.gz files."}
    }
}