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
        organism = "String"
        genome_fasta = "File"
        genes_gtf = "File"
    }

    String genome_name = if organism == "human" then "hg38" else "mm10"

    Int threads = 16
    Int mem_gb = 48
    Int disk_size = 50

    command <<<
        set -euo pipefail

        mkdir ref_dir
        
        split-pipe \
            --mode mkref \
            --genome_name ~{genome_name} \
            --fasta ~{genome_fasta} \
            --genes ~{genes_gtf} \
            --output_dir ref_dir

        tar -czf ref_dir.tar.gz ref_dir
    >>>

    output {
        File ref_dir_tar_gz = "ref_dir.tar.gz"

    }

    runtime {
        docker: "~{container_registry}/split-pipe:1.8.2"
        cpu: threads
        # cpuPlatform: "Intel Cascade Lake"  # necessary?
        memory: "{~mem_gb} GB"
        disks "local-disk ~{disk_size} HDD"
        bootDiskSizeGb: 40
        zones: zones

    }

    meta {
        description = "Make reference genome index from genome FASTA and GTF files."
    }

    parameter_meta {
        organism: {help: "Organism; used to select which reference genome is used. Either 'human' (HG38) or 'mouse' (MM10)."}
        genome_fasta: {help: "Path to the reference genome FASTA file."}
        genes_gtf: {help: "Path to the reference gene annotation GTF file."}
    }

}