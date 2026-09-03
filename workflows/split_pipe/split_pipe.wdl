version 1.0

# Call split-pipe to generate a matrix of counts from paired-end FASTQ files

import "../structs.wdl"


workflow split_pipe_wf {
    input {

    }

    output {
        
    }

    meta {
        description = "TODO!"
    }

    parameter_meta {
    
    }

}


task make_reference {
    input {
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
        # cpuPlatform: "Intel Cascade Lake"  # TODO: necessary?
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} HDD"
        bootDiskSizeGb: 40
        zones: zones
    }

    meta {
        description = "Make reference genome index from FASTA and GTF files."
    }

    parameter_meta {
        organism: {help: "Organism; used to select which reference genome is used. Either 'human' (HG38) or 'mouse' (MM10)."}
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
        Array[File] fastqs_R1s
        Array[File] fastqs_R2s
        File sample_table
        File? parameter_file
    }

    command <<<
        set -euo pipefail

        # Extract the genome reference dir and make the output sublibrary dir
        tar -xzf ~{genome_ref_tar_gz}
        mkdir ~{sublibrary_id}

        split-pipe \
            --mode ~{mode} \
            --chemistry ~{chemistry} \
            --kit ~{kit} \
            --fq1 ~{fastqs_R1s} \
            --fq2 ~{fastqs_R2s} \
            --genome_dir genome_ref \
            --output_dir ~{sublibrary_id} \
            --samp_list ~{sample_table} \
            --parfile ~{parameter_file}                 
            --dryrun      # To be removed after testing

        tar -czf ~{sublibrary_id}.tar.gz ~{sublibrary_id}        
    >>>

    output {
        File sublibrary_tar_gz = "~{sublibrary_id}.tar.gz"
    }

    runtime {
        docker: "~{container_registry}/split-pipe:1.8.2"
        cpu: threads
        # cpuPlatform: "Intel Cascade Lake"  # TODO: necessary?
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
    }
    # File of sublib dir paths to be fed as input to split-pipe combine
    File sublib_list = write_lines(sublibrary_ids)

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

    tar -czf combined_sublibs.tar.gz combined_sublibs
    >>>

    output {
        File combined_sublibs_tar_gz = "combined_sublibs.tar.gz"
    }

    meta {
        description = "Combine multiple sublibraries into a one combined object."
    }








     