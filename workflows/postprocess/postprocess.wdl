version 1.0
# Wrapper WDL for testing the Vireo and QC steps by pointing to intermediate combined sublibraries

import "../structs.wdl"

workflow postprocess_entrypoint {
	input {
		File combined_sublibs_tar_gz
        String dataset_id
        String team_id
        String kit
        String batch
        Array[ParseSample] parse_samples
        File sample_csv
        File subject_csv
        File vireo_join_script  # This will be baked into the image once verified
        File prep_metadata_script   # This will be baked into the image once verified
        String container_registry
        String zones = "us-central1-c us-central1-f"
	}
    
     scatter (parse_sample in parse_samples) {
        call prep_parse_sample {
            input:
                combined_sublibs_tar_gz = combined_sublibs_tar_gz,
                dataset_id = dataset_id,
                team_id = team_id,
                kit = kit,
                batch = batch,
                sample_name = parse_sample.parse_sample_name,
                donor_id_vireo_mapping = parse_sample.vireo_donor_ids,
                subject_id = parse_sample.subject_id,
                sample_csv = sample_csv,
                subject_csv = subject_csv,
                vireo_join_script = vireo_join_script,
                prep_metadata_script = prep_metadata_script,
                container_registry = container_registry,
                zones = zones
        }
    }
    
    output {
        Array[File] vireo_joined_cell_meta = prep_parse_sample.vireo_joined_cell_meta
        Array[File] prepped_adata = prep_parse_sample.prepped_adata
    }
}


task prep_parse_sample {
    input {
        File combined_sublibs_tar_gz
        String dataset_id
        String team_id
        String kit
        String batch
        String sample_name
        File? donor_id_vireo_mapping
        String? subject_id
        File sample_csv
        File subject_csv
        File vireo_join_script
        File prep_metadata_script
        String container_registry
        String zones
    }

    String prefix = "~{dataset_id}.~{kit}.~{sample_name}"
    String joined_meta = "~{prefix}_vireo_joined_cell_meta.csv"
    String adata_out = "~{prefix}_prepped.h5ad"

    Int threads = 4
    Int mem_gb = 16  
    Int disk_size = ceil((size(combined_sublibs_tar_gz, "GB") * 3) + 20)

    command <<<
        set -euo pipefail

        tar -xzf ~{combined_sublibs_tar_gz}
        dge_dir="combined_sublibs/~{sample_name}/DGE_filtered"

        python3 ~{vireo_join_script} \
            --cell_meta "${dge_dir}/cell_metadata.csv" \
            ~{"--vireo " + donor_id_vireo_mapping} \
            ~{"--subject_id " + subject_id} \
            --output ~{joined_meta}

        python3 ~{prep_metadata_script} \
            --data_dir "${dge_dir}" \
            --vireo_metadata_path ~{joined_meta} \
            --sample_path ~{sample_csv} \
            --subject_path ~{subject_csv} \
            --batch ~{batch} \
            --team ~{team_id} \
            --dataset ~{dataset_id} \
            --output ~{adata_out}
    >>>

    output {
        File vireo_joined_cell_meta = joined_meta  # TODO: necessary to keep?
        File prepped_adata = adata_out
    }

    runtime {
        docker: "~{container_registry}/sc_tools:1.3.0"
        cpu: threads
        cpuPlatform: "Intel Cascade Lake"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} SDD"
        preemptible: 3
        zones: zones
    }

    meta {
        description: "Per Parse sample: assign donors (vireo table or single known donor), compute QC metrics and scrublet scores, and attach ASAP sample/subject metadata."
    }

    parameter_meta {
        combined_sublibs_tar_gz: {help: "Tarball of the split-pipe comb output directory."}
		dataset_id: {help: "ASAP dataset identifier; used in output names and batch_id."}
		team_id: {help: "ASAP team identifier; used in batch_id."}
		kit: {help: "Parse kit run (e.g. PB-001); used in output names."}
		batch: {help: "SAMPLE.csv batch value for this kit run (e.g. Batch1)."}
		sample_name: {help: "Parse sample name; must match a comb output subdirectory."}
		donor_id_vireo_mapping: {help: "Vireo donor assignment TSV for a pooled Parse sample. Mutually exclusive with subject_id."}
		subject_id: {help: "Single known donor for an un-pooled Parse sample. Mutually exclusive with donor_id_vireo_mapping."}
		sample_csv: {help: "ASAP SAMPLE.csv."}
		subject_csv: {help: "ASAP SUBJECT.csv."}
	}
}