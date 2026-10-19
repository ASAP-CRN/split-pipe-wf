version 1.0

# Take the combined sublibraries output from split-pipe and build per-run filtered AnnData objects.

import "../structs.wdl"

workflow postprocess {
    input {
        File combined_sublibs_tar_gz
        String dataset_id
        String team_id
        String run_id
        String batch
        Array[ParseSample] parse_samples
        File sample_csv
        File subject_csv

        # Filtering parameters
        Int pct_counts_mt_max
        Float doublet_score_max
        Array[Int] total_counts_limits
        Array[Int] n_genes_by_counts_limits
        Float n_mads_lower
        Float n_mads_upper

        String container_registry
        String zones = "us-central1-c us-central1-f"
    }

    scatter (parse_sample in parse_samples) {
        call prep_parse_sample {
            input:
                combined_sublibs_tar_gz = combined_sublibs_tar_gz,
                dataset_id = dataset_id,
                team_id = team_id,
                run_id = run_id,
                batch = batch,
                sample_name = parse_sample.parse_sample_name,
                donor_id_vireo_mapping = parse_sample.vireo_donor_ids,
                subject_id = parse_sample.subject_id,
                sample_csv = sample_csv,
                subject_csv = subject_csv,
                container_registry = container_registry,
                zones = zones
        }
    }

    call merge_and_plot_qc_metrics {
        input:
            dataset_id = dataset_id,
            run_id = run_id,
            prepped_adata_objects = prep_parse_sample.prepped_adata,
            container_registry = container_registry,
            zones = zones
    }

    call filter {
        input:
            dataset_id = dataset_id,
            run_id = run_id,
            merged_adata_object = merge_and_plot_qc_metrics.merged_adata_object,
            pct_counts_mt_max = pct_counts_mt_max,
            doublet_score_max = doublet_score_max,
            total_counts_limits = total_counts_limits,
            n_genes_by_counts_limits = n_genes_by_counts_limits,
            n_mads_lower = n_mads_lower,
            n_mads_upper = n_mads_upper,
            container_registry = container_registry,
            zones = zones
    }

    output {
        Array[File] prepped_adata = prep_parse_sample.prepped_adata
        File merged_adata = merge_and_plot_qc_metrics.merged_adata_object
        File initial_metadata = merge_and_plot_qc_metrics.qc_initial_metadata_csv
        Array[File] qc_plots = merge_and_plot_qc_metrics.qc_plots_png
        File filtered_adata_object = filter.filtered_adata_object
        File unassigned_removed_cells_csv = filter.unassigned_removed_cells_csv
        File filtered_out_cells_csv = filter.filtered_out_cells_csv
    }

    parameter_meta {
        combined_sublibs_tar_gz: {help: "Tarball of the split-pipe comb output directory for one Parse run."}
        dataset_id: {help: "ASAP dataset identifier; used in output names and batch_id."}
        team_id: {help: "ASAP team identifier; used in batch_id."}
        run_id: {help: "Parse run ID (e.g. PB-001); used in output names and cell-ID suffixes."}
        batch: {help: "SAMPLE.csv batch value for this Parse run (e.g. Batch1)."}
        parse_samples: {help: "Parse samples to process; each has a vireo donor table (pooled) or a single subject_id (un-pooled)."}
        sample_csv: {help: "ASAP SAMPLE.csv."}
        subject_csv: {help: "ASAP SUBJECT.csv."}
        pct_counts_mt_max: {help: "Maximum percentage of mitochondrial gene counts allowed per cell. [10]"}
        doublet_score_max: {help: "Maximum doublet detection score threshold. [0.2]"}
        total_counts_limits: {help: "Absolute minimum and maximum total UMI counts per cell, applied on top of the MAD-based thresholds. [500, 150000]"}
        n_genes_by_counts_limits: {help: "Absolute minimum and maximum number of genes detected per cell, applied on top of the MAD-based thresholds. [300, 15000]"}
        n_mads_lower: {help: "Number of MADs below the per-sample median allowed for total counts and genes detected. [3]"}
        n_mads_upper: {help: "Number of MADs above the per-sample median allowed for total counts and genes detected. [5]"}
        container_registry: {help: "Container registry where workflow Docker images are hosted."}
        zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
    }
}


task prep_parse_sample {
    input {
        File combined_sublibs_tar_gz
        String dataset_id
        String team_id
        String run_id
        String batch
        String sample_name
        File? donor_id_vireo_mapping
        String? subject_id
        File sample_csv
        File subject_csv
        String container_registry
        String zones
    }

    String prefix = "~{dataset_id}.~{run_id}.~{sample_name}"
    String donor_assignments = "~{prefix}.donor_assignments.csv"
    String adata_out = "~{prefix}.prepped.h5ad"

    Int threads = 4
    Int mem_gb = 16
    Int disk_size = ceil((size(combined_sublibs_tar_gz, "GB") * 3) + 20)

    command <<<
        set -euo pipefail

        tar -xzf ~{combined_sublibs_tar_gz}
        dge_dir="combined_sublibs/~{sample_name}/DGE_filtered"

        assign_donors.py \
            --cell-meta "${dge_dir}/cell_metadata.csv" \
            ~{"--vireo " + donor_id_vireo_mapping} \
            ~{"--subject-id " + subject_id} \
            --output ~{donor_assignments}

        prep_parse_sample.py \
            --data-dir "${dge_dir}" \
            --donor-assignments ~{donor_assignments} \
            --sample-path ~{sample_csv} \
            --subject-path ~{subject_csv} \
            --batch ~{batch} \
            --team ~{team_id} \
            --dataset ~{dataset_id} \
            --output ~{adata_out}
    >>>

    output {
        File vireo_joined_cell_meta = donor_assignments
        File prepped_adata = adata_out
    }

    runtime {
        docker: "~{container_registry}/parsebio_tools:0.1.0"
        cpu: threads
        cpuPlatform: "Intel Cascade Lake"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} SSD"
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
        run_id: {help: "Parse run ID (e.g. PB-001); used in output names."}
        batch: {help: "SAMPLE.csv batch value for this Parse run (e.g. Batch1)."}
        sample_name: {help: "Parse sample name; must match a comb output subdirectory."}
        donor_id_vireo_mapping: {help: "Vireo donor assignment TSV for a pooled Parse sample. Mutually exclusive with subject_id."}
        subject_id: {help: "Single known donor for an un-pooled Parse sample. Mutually exclusive with donor_id_vireo_mapping."}
        sample_csv: {help: "ASAP SAMPLE.csv."}
        subject_csv: {help: "ASAP SUBJECT.csv."}
        container_registry: {help: "Container registry where workflow Docker images are hosted."}
        zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
    }
}


task merge_and_plot_qc_metrics {
    input {
        String dataset_id
        String run_id
        Array[File] prepped_adata_objects
        String container_registry
        String zones
    }

    String prefix = "~{dataset_id}.~{run_id}"

    Int mem_gb = ceil(size(prepped_adata_objects, "GB") * 2.4 + 20)
    Int disk_size = ceil(size(prepped_adata_objects, "GB") * 3 + 50)

    command <<<
        set -euo pipefail

        while read -r adata_object || [[ -n "${adata_object}" ]]; do
            adata_path=$(realpath "${adata_object}")
            sample=$(basename "${adata_path}" ".prepped.h5ad")
            sample="${sample#~{dataset_id}.}"
            echo -e "${sample}\t${adata_path}" >> adata_samples_paths.tsv
        done < ~{write_lines(prepped_adata_objects)}

        merge_and_plot_qc \
            --adata-objects-fofn adata_samples_paths.tsv \
            --adata-output ~{prefix}.merged_unfiltered.h5ad \
            --output-metadata-file ~{prefix}.initial_metadata.csv

        mv "plots/violin_n_genes_by_counts.png" "plots/~{prefix}.n_genes_by_counts.violin.png"
		mv "plots/violin_total_counts.png" "plots/~{prefix}.total_counts.violin.png"
		mv "plots/violin_pct_counts_mt.png" "plots/~{prefix}.pct_counts_mt.violin.png"
		mv "plots/violin_pct_counts_rb.png" "plots/~{prefix}.pct_counts_rb.violin.png"
		mv "plots/violin_doublet_score.png" "plots/~{prefix}.doublet_score.violin.png"
        
    >>>

    output {
        File merged_adata_object = "~{prefix}.merged_unfiltered.h5ad"
        File qc_initial_metadata_csv = "~{prefix}.initial_metadata.csv"
        Array[File] qc_plots_png = [
            "plots/~{prefix}.n_genes_by_counts.violin.png",
            "plots/~{prefix}.total_counts.violin.png",
            "plots/~{prefix}.pct_counts_mt.violin.png",
            "plots/~{prefix}.pct_counts_rb.violin.png",
            "plots/~{prefix}.doublet_score.violin.png"
        ]
    }

    runtime {
        docker: "~{container_registry}/parsebio_tools:0.1.0"
        cpu: 4
        cpuPlatform: "Intel Cascade Lake"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 3
        bootDiskSizeGb: 40
        zones: zones
    }

    meta {
        description: "Merges Parse Bio sample-level AnnData objects to a single AnnData object and generates pre-filtering QC violin plots for key metrics."
    }

    parameter_meta {
        dataset_id: {help: "ASAP dataset identifier; used in output names and batch_id."}
        run_id: {help: "Parse run ID (e.g. PB-001); used in output names."}
        prepped_adata_objects: {help: "An array of Parse Sample-level AnnData objects to merge."}
        container_registry: {help: "Container registry where workflow Docker images are hosted."}
        zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
    }
}


task filter {
    input {
        String dataset_id
        String run_id
        File merged_adata_object

        Int pct_counts_mt_max
        Float doublet_score_max
        Array[Int] total_counts_limits
        Array[Int] n_genes_by_counts_limits
        Float n_mads_lower
        Float n_mads_upper

        String container_registry
        String zones
    }

    String prefix = "~{dataset_id}.~{run_id}"

    Int calc_mem_gb = ceil(size(merged_adata_object, "GB") * 18 + 20)
    Int mem_gb = if calc_mem_gb > 624 then 624 else calc_mem_gb
    Int disk_size = ceil(size(merged_adata_object, "GB") * 4 + 20)

    command <<<
        set -euo pipefail

        drop_unassigned_cells.py \
            --adata-input ~{merged_adata_object} \
            --adata-output ~{prefix}.assigned.h5ad \
            --output-removed-cells-file ~{prefix}.unassigned_removed_cells.csv

        filter \
            --adata-input ~{prefix}.assigned.h5ad \
            --pct-counts-mt-max ~{pct_counts_mt_max} \
            --doublet-score-max ~{doublet_score_max} \
            --total-counts-limits ~{sep=' ' total_counts_limits} \
            --n-genes-by-counts-limits ~{sep=' ' n_genes_by_counts_limits} \
            --n-mads-lower ~{n_mads_lower} \
            --n-mads-upper ~{n_mads_upper} \
            --adata-output ~{prefix}.filtered.h5ad \
            --output-removed-cells-file ~{prefix}.filtered_out_cells.csv
    >>>

    output {
        File filtered_adata_object = "~{prefix}.filtered.h5ad"
        File unassigned_removed_cells_csv = "~{prefix}.unassigned_removed_cells.csv"
        File filtered_out_cells_csv = "~{prefix}.filtered_out_cells.csv"
    }

    runtime {
        docker: "~{container_registry}/parsebio_tools:0.1.0"
        cpu: 4
        cpuPlatform: "Intel Cascade Lake"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 3
        bootDiskSizeGb: 40
        zones: zones
    }

    meta {
        description: "Filters low-quality cells from the merged AnnData object based on mitochondrial content, doublet score, total UMI counts, and number of detected genes. UMI count and detected gene thresholds combine absolute limits with per-sample median absolute deviation (MAD) cutoffs. Lists every removed cell and the filters it failed."
    }

    parameter_meta {
        dataset_id: {help: "ASAP dataset identifier; used in output names and batch_id."}
        run_id: {help: "Parse run ID (e.g. PB-001); used in output names."}
        merged_adata_object: {help: "Merged AnnData object."}
        pct_counts_mt_max: {help: "Maximum percentage of mitochondrial gene counts allowed per cell. [10]"}
        doublet_score_max: {help: "Maximum doublet detection score threshold. [0.2]"}
        total_counts_limits: {help: "Absolute minimum and maximum total UMI (unique molecular identifier) counts per cell; applied on top of the MAD-based thresholds. [500, 150000]"}
        n_genes_by_counts_limits: {help: "Absolute minimum and maximum number of genes detected per cell (genes with at least one count); applied on top of the MAD-based thresholds. [300, 15000]"}
        n_mads_lower: {help: "Number of median absolute deviations below the per-sample median allowed for total UMI counts and number of genes detected per cell. [3]"}
        n_mads_upper: {help: "Number of median absolute deviations above the per-sample median allowed for total UMI counts and number of genes detected per cell. [5]"}
        container_registry: {help: "Container registry where workflow Docker images are hosted."}
        zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
    }
}