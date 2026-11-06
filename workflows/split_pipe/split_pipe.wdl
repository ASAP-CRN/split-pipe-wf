version 1.0

# Call split-pipe to generate a matrix of counts from paired-end FASTQ files

import "../structs.wdl"

workflow split_pipe {
    input {
        String asap_dataset_id
        String run_id
        Array[Sublibrary] sublibraries
        String genome_name
        File? genome_ref_tar_gz
        File genome_fasta
        File genes_gtf

        String chemistry
        String kit
        File sample_list_tsv
        File? parameter_file
        
        String workflow_name
        String workflow_version
        String workflow_release
        String run_timestamp
        
        String raw_data_path_prefix
        String billing_project
        String container_registry
        String zones
    }

    String sub_workflow_name = "split_pipe"
    String sub_workflow_version = "1.0.0"
    String raw_data_path = "~{raw_data_path_prefix}/~{sub_workflow_name}/~{sub_workflow_version}/~{run_timestamp}"

    Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]
    
    # Build reference genome index if it was not provided
    if (!defined(genome_ref_tar_gz)) {
        call make_reference {
        input:
            genome_name = genome_name,
            genome_fasta = genome_fasta,
            genes_gtf = genes_gtf,
            container_registry = container_registry,
            zones = zones
        }
    }

    File genome_ref_final = select_first([
        genome_ref_tar_gz, 
        make_reference.genome_ref_tar_gz
    ])

    # Run split-pipe on each sublibrary
    scatter (sublibrary in sublibraries) {
        String sublibrary_id = sublibrary.sublibrary_id  # Declare so it can be gathered for array input for combine_sublibraries
        
        call split_pipe_sublibrary {
            input:
                genome_ref_tar_gz = genome_ref_final,
                asap_dataset_id = asap_dataset_id,
                sublibrary_id = sublibrary_id,
                chemistry = chemistry,
                kit = kit,
                fastq_R1s = sublibrary.fastq_R1s,
                fastq_R2s = sublibrary.fastq_R2s,
                sample_list_tsv = sample_list_tsv,
                parameter_file = parameter_file,
                raw_data_path = raw_data_path,
                workflow_info = workflow_info,
                billing_project = billing_project,
                container_registry = container_registry,
                zones = zones
        }
    }

    Array[String] sublibrary_ids = sublibrary_id

    # Combine the output sublibaries into one object
    call combine_sublibraries {
        input:
            asap_dataset_id = asap_dataset_id,
            run_id = run_id,
            genome_ref_tar_gz = genome_ref_final,
            sublibrary_tars = split_pipe_sublibrary.sublibrary_tar_gz,  # Gathered Array[File] from the split-pipe scatter
            sublibrary_ids = sublibrary_ids,
            parameter_file = parameter_file,
            raw_data_path = raw_data_path,
            workflow_info = workflow_info,
            billing_project = billing_project,
            container_registry = container_registry,
            zones = zones
    }

    output {
        Array[File] sublibrary_tar_gz = split_pipe_sublibrary.sublibrary_tar_gz  # Intermediate, not uploaded
        Array[File] sublibrary_bams = split_pipe_sublibrary.bam #!FileCoercion
        Array[File] sublibrary_star_logs = split_pipe_sublibrary.star_log #!FileCoercion
        File combined_sublibs_tar_gz = combine_sublibraries.combined_sublibs_tar_gz #!FileCoercion
        File all_sample_summary_html = combine_sublibraries.all_sample_summary_html #!FileCoercion
        File agg_sample_summary_csv = combine_sublibraries.agg_sample_summary_csv #!FileCoercion
    }

    meta {
        description: "This workflow creates the indexed reference genome if it is not provided, runs split-pipe on each sublibrary and combines the results into a single object."
    }

    parameter_meta {
        asap_dataset_id: {help: "ASAP dataset ID (e.g. DS_PMDBS_0038); used in output names."}
        run_id: {help: "Parse run ID (e.g. PB-001); used in output names and cell-ID suffixes."}
        sublibraries: {help: "Array of sublibrary structs containing sublibrary IDs and FASTQ file paths."}
        genome_name: {help: "The name of which reference genome is used. Either 'hg38' (human) or 'mm10' (mouse)."}
        genome_ref_tar_gz: {help: "Optional prebuilt split-pipe mkref index (tar.gz). If omitted, one is built from genome_fasta and genes_gtf."}
        genome_fasta: {help: "Path to the reference genome FASTA file."}
        genes_gtf: {help: "Path to the reference gene annotation GTF file."}
        chemistry: {help: "Parse Bio chemistry: v1-v4"}
        kit: {help: "Parse Bio kit: WT_mini, WT, WT_mega, WT_mega_384, WT_penta, WT_penta_384"}
        sample_list_tsv: {help: "Parse sample list (--samp_list): one line per sample, '<sample_name> <round-1 wells>', e.g. 'pool_1 A1-B3'. Shared by every sublibrary in the run."}
        parameter_file: {help: "Optional split-pipe --parfile ('key value' per line), passed to every sublibrary and comb call. Use split_pipe_params.txt (log_memory_use False) to avoid a split-pipe memory-profiler race that can fail a run."}
        raw_data_path_prefix: {help: "Raw data bucket prefix (`<raw_data_bucket>/workflow_execution/<workflow_name>`); outputs are uploaded to `<prefix>/split_pipe/<sub_workflow_version>/<run_timestamp>`."}
        workflow_name: {help: "Workflow name; stored in the file-level manifest and final manifest with all saved files."}
        workflow_version: {help: "Workflow version; stored in the file-level manifest and final manifest with all saved files."}
        workflow_release: {help: "GitHub release; stored in the file-level manifest and final manifest with all saved files."}
        run_timestamp: {help: "UTC timestamp; stored in the file-level manifest and final manifest with all saved files."}
        billing_project: {help: "Billing project to charge GCP costs."}
        container_registry: {help: "Container registry where workflow Docker images are hosted."}
        zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
    }
}


task make_reference {
    input {
        String genome_name
        File genome_fasta
        File genes_gtf
        String container_registry
        String zones
    }

    String annotation_label = basename(genes_gtf, ".genes.gtf.gz")
    String genome_ref_final = "split_pipe_1.8.2.~{genome_name}.~{annotation_label}.genome_ref.tar.gz"

    Int threads = 16
    Int mem_gb = 48
    Int disk_size = 400

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

        tar -czf "~{genome_ref_final}" genome_ref
    >>>

    output {
        File genome_ref_tar_gz = genome_ref_final
    }

    runtime {
        docker: "~{container_registry}/split_pipe:1.8.2"
        cpu: threads
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 3
        zones: zones
    }

    meta {
        description: "Make reference genome index from FASTA and GTF files."
    }

    parameter_meta {
        genome_name: {help: "The name of which reference genome is used. Either 'hg38' (human) or 'mm10' (mouse)."}
        genome_fasta: {help: "Path to the reference genome FASTA file."}
        genes_gtf: {help: "Path to the reference gene annotation GTF file."}
        container_registry: {help: "Container registry where workflow Docker images are hosted."}
        zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
    }
}


task split_pipe_sublibrary {
    input {
        String asap_dataset_id
        String sublibrary_id
        File genome_ref_tar_gz
        String chemistry
        String kit
        Array[File] fastq_R1s
        Array[File] fastq_R2s
        File sample_list_tsv
        File? parameter_file
        
        String raw_data_path
        Array[Array[String]] workflow_info
        String billing_project
        String container_registry
        String zones
    }

    Int threads = 16
    Int mem_gb = 32
    Int disk_size = ceil((size(genome_ref_tar_gz, "GB") + (size(flatten([fastq_R1s, fastq_R2s]), "GB")) * 5) + 50 )

    # Per sublibrary-specific deliverables - the rest is carried forward in combine_sublibraries
    String bam_name = "~{asap_dataset_id}.~{sublibrary_id}.name_sorted.bam"
    String star_log_name = "~{asap_dataset_id}.~{sublibrary_id}.barcode_headLog.final.out"

    command <<<
        set -euo pipefail

        # Extract the genome reference dir and make the output sublibrary dir
        mkdir genome_ref
        tar -xzf ~{genome_ref_tar_gz} -C genome_ref --strip-components=1
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
            ~{if defined(parameter_file) then "--parfile " + select_first([parameter_file]) else ""}

        # split-pipe can exit 0 having produced nothing; require DGE output for every listed sample
        while read -r sample _; do
            [[ -z "${sample}" ]] && continue
            if [[ ! -s "~{sublibrary_id}/${sample}/DGE_filtered/count_matrix.mtx" ]]; then
                echo "ERROR: split-pipe produced no DGE output for sample ${sample}" >&2
                exit 1
            fi
        done < ~{sample_list_tsv}

        # Name-sorted, tag-annotated BAM (not needed by comb); moved out so it is not carried in the tarball
        mv ~{sublibrary_id}/process/barcode_headAligned_anno.bam "~{bam_name}"
        cp ~{sublibrary_id}/process/barcode_headLog.final.out "~{star_log_name}"

        upload_outputs \
            -b ~{billing_project} \
            -d ~{raw_data_path} \
            -i ~{write_tsv(workflow_info)} \
            -o "~{bam_name}" \
            -o "~{star_log_name}"

        tar -czf ~{asap_dataset_id}.~{sublibrary_id}.tar.gz ~{sublibrary_id}
    >>>

    output {
        File sublibrary_tar_gz = "~{asap_dataset_id}.~{sublibrary_id}.tar.gz"  # Passed to comb, not uploaded b/c it becomes part of comb output
        String bam = "~{raw_data_path}/~{bam_name}"
        String star_log = "~{raw_data_path}/~{star_log_name}"
    }

    runtime {
        docker: "~{container_registry}/split_pipe:1.8.2"
        cpu: threads
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} SSD"
        zones: zones
    }

    meta {
        description: "Run split-pipe on mode: 'all' for the given sublibrary."
    }

    parameter_meta {
        asap_dataset_id: {help: "ASAP dataset ID (e.g. DS_PMDBS_0038); used in output names."}
        sublibrary_id: {help: "Identifier for the sublibrary being processed."}
        genome_ref_tar_gz: {help: "split-pipe mkref index (tar.gz), from make_reference or supplied to the workflow."}
        chemistry: {help: "Parse Bio chemistry: v1-v4"}
        kit: {help: "Parse Bio kit: WT_mini, WT, WT_mega, WT_mega_384, WT_penta, WT_penta_384"}
        fastq_R1s: {help: "Array of FASTQ files for read 1."}
        fastq_R2s: {help: "Array of FASTQ files for read 2."}
        sample_list_tsv: {help: "Parse sample list (--samp_list): one line per sample, '<sample_name> <round-1 wells>', e.g. 'pool_1 A1-B3'. Shared by every sublibrary in the run."}
        parameter_file: {help: "Optional split-pipe --parfile; 'key value' per line."}
        raw_data_path: {help: "Raw data bucket path to upload outputs to."}
        workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
        billing_project: {help: "Billing project to charge GCP costs."}
        container_registry: {help: "Container registry where workflow Docker images are hosted."}
        zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
    }
}


task combine_sublibraries {
    input {
        String asap_dataset_id
        String run_id
        File genome_ref_tar_gz
        Array[File] sublibrary_tars
        Array[String] sublibrary_ids
        File? parameter_file

        String raw_data_path
        Array[Array[String]] workflow_info
        String billing_project
        String container_registry
        String zones
    }
    
    File sublib_list = write_lines(sublibrary_ids)
    String comb_name = "~{asap_dataset_id}.~{run_id}.combined_sublibs"
    String comb_tar_gz = "~{comb_name}.tar.gz"
    
    Int threads = 16
    Int mem_gb = 96
    Int disk_size = ceil((size(sublibrary_tars, "GB") * 8) + (size(genome_ref_tar_gz, "GB") * 3) + 50)

    command <<<
        set -euo pipefail
        
        # comb needs the genome ref around, but calling it as an argument leads to a tricky 
        # interaction with cromwell pathing. Having it sit in the CWD gets around this.
        mkdir genome_ref
        tar -xzf ~{genome_ref_tar_gz} -C genome_ref --strip-components=1

        # Extract each sublibrary .tar.gz into sublibs/
        mkdir -p sublibs
        for sublib_tar in ~{sep=" " sublibrary_tars}; do
            tar -xzf "${sublib_tar}" -C sublibs
        done

        /usr/bin/time -v \
        split-pipe \
            --mode comb \
            --sublib_list ~{sublib_list} \
            --sublib_pref "$PWD/sublibs/" \
            --output_dir combined_sublibs \
            --nthreads ~{threads} \
            ~{if defined(parameter_file) then "--parfile " + select_first([parameter_file]) else ""}

        # split-pipe can exit 0 having produced nothing
        if [[ ! -s combined_sublibs/all-sample/DGE_filtered/count_matrix.mtx ]]; then
            echo "ERROR: split-pipe comb produced no DGE output" >&2
            exit 1
        fi

        # Copy the summary reports to view alongside the compressed output
        cp combined_sublibs/all-sample_analysis_summary.html "~{comb_name}.all-sample_analysis_summary.html"
        cp combined_sublibs/agg_sample_summary.csv "~{comb_name}.agg_sample_summary.csv"
        tar -czf "~{comb_tar_gz}" combined_sublibs

        upload_outputs \
            -b ~{billing_project} \
            -d ~{raw_data_path} \
            -i ~{write_tsv(workflow_info)} \
            -o "~{comb_name}.all-sample_analysis_summary.html" \
            -o "~{comb_name}.agg_sample_summary.csv" \
            -o "~{comb_tar_gz}"
    >>>

    output {
        String combined_sublibs_tar_gz = "~{raw_data_path}/~{comb_tar_gz}"
        String all_sample_summary_html = "~{raw_data_path}/~{comb_name}.all-sample_analysis_summary.html"
        String agg_sample_summary_csv = "~{raw_data_path}/~{comb_name}.agg_sample_summary.csv"
    }

    runtime {
        docker: "~{container_registry}/split_pipe:1.8.2"
        cpu: threads
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} SSD"
        zones: zones
    }

    meta {
        description: "Combine multiple sublibraries into a single combined object."
    }

    parameter_meta {
        asap_dataset_id: {help: "ASAP dataset ID (e.g. DS_PMDBS_0038); used in output names."}
        run_id: {help: "Parse run ID (e.g. PB-001); used in output names and cell-ID suffixes."}
        genome_ref_tar_gz: {help: "split-pipe mkref index (tar.gz), from make_reference or supplied to the workflow."}
        sublibrary_tars: {help: "Array of sublibrary tar.gz files to be combined."}
        sublibrary_ids: {help: "Array of sublibrary identifiers corresponding to the tar.gz files."}
        parameter_file: {help: "Optional split-pipe --parfile; 'key value' per line."}
        raw_data_path: {help: "Raw data bucket path to upload outputs to."}
        workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
        billing_project: {help: "Billing project to charge GCP costs."}
        container_registry: {help: "Container registry where workflow Docker images are hosted."}
        zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
    }
}