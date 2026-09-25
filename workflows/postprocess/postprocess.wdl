version 1.0
# Wrapper WDL for testing the Vireo and QC steps by pointing to intermediate combined sublibraries

import "../structs.wdl"

workflow preprocess_entrypoint {
	input {
		File combined_sublibs_tar_gz
        String dataset_id
        Array[ParseSample] parse_samples
        File? donor_id_vireo_mapping
        File vireo_join_script  # This will be baked into the image once verifiedv
        String container_registry
        String zones = "us-central1-c us-central1-f"
	}
    
     scatter (parse_sample in parse_samples) {
        call vireo_and_qc_metrics {
            input:
                combined_sublibs_tar_gz = combined_sublibs_tar_gz,
                dataset_id = dataset_id,
                sample_name = parse_sample.parse_sample_name,
                donor_id_vireo_mapping = donor_id_vireo_mapping,
                vireo_join_script = vireo_join_script,
                container_registry = container_registry,
                zones = zones
        }
    }
    
    output {
        Array[File] vireo_joined_cell_meta = vireo_and_qc_metrics.vireo_joined_cell_meta
    }
}


task vireo_and_qc_metrics {
    input {
        File combined_sublibs_tar_gz
        String dataset_id
        String sample_name
        File? donor_id_vireo_mapping
        File vireo_join_script  # This will be baked into the image once verified
        String container_registry
        String zones
    }

    Int threads = 4
    Int mem_gb = 8  
    Int disk_size = ceil((size(combined_sublibs_tar_gz, "GB") * 3) + 10)

    command <<<
        set -euo pipefail

        mkdir combined_sublibs
        tar -xzf ~{combined_sublibs_tar_gz} -C combined_sublibs # TODO: dataset_id name
        cell_meta="combined_sublibs/~{sample_name}/DGE_unfiltered/cell_metadata.csv"

        if [ ! -f "$cell_meta" ]; then
            echo "Cell metadata file for sample $sample_name not found: $cell_meta"
            exit 1
        fi

        python3 ~{vireo_join_script} "$cell_meta" ~{donor_id_vireo_mapping} ~{dataset_id}.~{sample_name}_vireo_joined_cell_meta.csv
    >>>

    output {
        File vireo_joined_cell_meta = "~{dataset_id}.~{sample_name}_vireo_joined_cell_meta.csv"
    }

    runtime {
        docker: "~{container_registry}/sc_tools:1.3.0"
        cpu: threads
        cpuPlatform: "Intel Cascade Lake"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} HDD"
        preemptible: 3
        zones: zones
    }

    meta {
        description: "TODO"
    }

    parameter_meta {
        combined_sublibs_tar_gz: {help: "Tar.gz file containing combined sublibraries."}
        dataset_id: {help: "The unique identifier for the dataset."}
        sample_name: {help: "TODO"}
        donor_id_vireo_mapping: {help: "Mapping file of donor IDs deconvoluted by Vireo."}
    }
}