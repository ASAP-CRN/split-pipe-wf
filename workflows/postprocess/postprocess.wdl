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
        File? sample_map

        # Filtering parameters
        Int pct_counts_mt_max
        Float doublet_score_max
        Array[Int] total_counts_limits
        Array[Int] n_genes_by_counts_limits
        Float n_mads_lower
        Float n_mads_upper

        # Allen Institute's Map My Cells
        String mmc_taxonomy
        File allen_brain_mmc_precomputed_stats_h5
        File? allen_brain_mmc_marker_genes_json

        String container_registry
        String zones = "us-central1-c us-central1-f"
    }

    String run_prefix = "~{dataset_id}.~{run_id}"

    scatter (parse_sample in parse_samples) {

        # sample_prefix labels a parse sample, merge_key (dataset_id excluded for
        # brevity) labels cells within a parse sample so they are unique upon merge
        String sample_prefix = "~{run_prefix}.~{parse_sample.parse_sample_name}"
        String merge_key = "~{run_id}.~{parse_sample.parse_sample_name}"

        call prep_parse_sample {
            input:
                combined_sublibs_tar_gz = combined_sublibs_tar_gz,
                output_prefix = sample_prefix,
                dataset_id = dataset_id,
                team_id = team_id,
                batch = batch,
                sample_name = parse_sample.parse_sample_name,
                donor_id_vireo_mapping = parse_sample.vireo_donor_ids,
                subject_id = parse_sample.subject_id,
                sample_csv = sample_csv,
                subject_csv = subject_csv,
                sample_map = sample_map,
                container_registry = container_registry,
                zones = zones
        }
    }

    call merge_and_plot_qc_metrics {
        input:
            output_prefix = run_prefix,
            merge_keys = merge_key,
            prepped_adata_objects = prep_parse_sample.prepped_adata,
            container_registry = container_registry,
            zones = zones
    }

    call filter {
        input:
            output_prefix = run_prefix,
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

    call map_cell_types {
        input:
            output_prefix = run_prefix,
            filtered_adata_object = filter.filtered_adata_object,
            mmc_taxonomy = mmc_taxonomy,
            allen_brain_mmc_precomputed_stats_h5 = allen_brain_mmc_precomputed_stats_h5,
            allen_brain_mmc_marker_genes_json = allen_brain_mmc_marker_genes_json,
            container_registry = container_registry,
            zones = zones
    }

    call add_mapped_cell_types {
        input:
            output_prefix = run_prefix,
            filtered_adata_object = filter.filtered_adata_object,
            mmc_results_csv = map_cell_types.mmc_results_csv,
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
        File mmc_extended_results_json = map_cell_types.mmc_extended_results_json
        File mmc_results_csv = map_cell_types.mmc_results_csv
        File mmc_log_txt = map_cell_types.mmc_log_txt
        File mmc_adata_object = add_mapped_cell_types.mmc_adata_object
        File mmc_results_parquet = add_mapped_cell_types.mmc_results_parquet
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
        sample_map: {help: "CSV of parse_sample, subject_id, sample_id; needed when a subject has more than one sample in the batch."}
        pct_counts_mt_max: {help: "Maximum percentage of mitochondrial gene counts allowed per cell. [10]"}
        doublet_score_max: {help: "Maximum doublet detection score threshold. [0.2]"}
        total_counts_limits: {help: "Absolute minimum and maximum total UMI counts per cell, applied on top of the MAD-based thresholds. [500, 150000]"}
        n_genes_by_counts_limits: {help: "Absolute minimum and maximum number of genes detected per cell, applied on top of the MAD-based thresholds. [300, 15000]"}
        n_mads_lower: {help: "Number of MADs below the per-sample median allowed for total counts and genes detected. [3]"}
        n_mads_upper: {help: "Number of MADs above the per-sample median allowed for total counts and genes detected. [5]"}
        mmc_taxonomy: {help: "Cell type taxonomy of the precomputed stats reference; appended to MMC output filenames. Must match allen_brain_mmc_precomputed_stats_h5."}
        allen_brain_mmc_precomputed_stats_h5: {help: "A precomputed statistics file from the Allen Brain Cell Atlas containing reference statistics (the average gene expression profile per cell type cluster and cell type taxonomy)."}
        allen_brain_mmc_marker_genes_json: {help: "Optional marker gene lookup (parent node -> marker genes), e.g. extracted from a previous run's extended_results.json. If provided, skips on-the-fly marker computation."}
        container_registry: {help: "Container registry where workflow Docker images are hosted."}
        zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
    }
}


task prep_parse_sample {
    input {
        String output_prefix
        File combined_sublibs_tar_gz
        String dataset_id
        String team_id
        String batch
        String sample_name
        File? donor_id_vireo_mapping
        String? subject_id
        File sample_csv
        File subject_csv
        File? sample_map
        String container_registry
        String zones
    }

    String donor_assignments = "~{output_prefix}.donor_assignments.csv"
    String adata_out = "~{output_prefix}.prepped.h5ad"

    Int threads = 4
    Int mem_gb = 16
    Int disk_size = ceil(size(combined_sublibs_tar_gz, "GB") + 20)

    command <<<
        set -euo pipefail

        # Only need the filtered DGE output for this Parse sample
        dge_dir="combined_sublibs/~{sample_name}/DGE_filtered"
        tar -xzf ~{combined_sublibs_tar_gz} "${dge_dir}"

        /usr/bin/time -v \
        assign_donors.py \
            --cell-meta "${dge_dir}/cell_metadata.csv" \
            ~{"--vireo " + donor_id_vireo_mapping} \
            ~{"--subject-id " + subject_id} \
            --dataset-id ~{dataset_id} \
            --output ~{donor_assignments}

        /usr/bin/time -v \
        prep_parse_sample.py \
            --data-dir "${dge_dir}" \
            --donor-assignments ~{donor_assignments} \
            --sample-path ~{sample_csv} \
            --subject-path ~{subject_csv} \
            --batch ~{batch} \
            --team ~{team_id} \
            --dataset ~{dataset_id} \
            ~{"--sample-map " + sample_map} \
            --output ~{adata_out}
    >>>

    output {
        File donor_assignments_csv = donor_assignments
        File prepped_adata = adata_out
    }

    runtime {
        docker: "~{container_registry}/parsebio_tools:0.3.1"
        cpu: threads
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} SSD"
        bootDiskSizeGb: 40
        preemptible: 3
        zones: zones
    }

    meta {
        description: "Per Parse sample: assign donors (vireo table or single known donor), compute QC metrics and scrublet scores, and attach ASAP sample/subject metadata."
    }

    parameter_meta {
        output_prefix: {help: "Prefix for output files of the form <dataset_id>.<run_id>.<parse_sample>."}
        combined_sublibs_tar_gz: {help: "Tarball of the split-pipe comb output directory."}
        dataset_id: {help: "ASAP dataset identifier; used in output names and batch_id."}
        team_id: {help: "ASAP team identifier; used in batch_id."}
        batch: {help: "SAMPLE.csv batch value for this Parse run (e.g. Batch1)."}
        sample_name: {help: "Parse sample name; must match a comb output subdirectory."}
        donor_id_vireo_mapping: {help: "Vireo donor assignment TSV for a pooled Parse sample. Mutually exclusive with subject_id."}
        subject_id: {help: "Single known donor for an un-pooled Parse sample. Mutually exclusive with donor_id_vireo_mapping."}
        sample_csv: {help: "ASAP SAMPLE.csv."}
        subject_csv: {help: "ASAP SUBJECT.csv."}
        sample_map: {help: "CSV of parse_sample, subject_id, sample_id; needed when a subject has more than one sample in the batch."}
        container_registry: {help: "Container registry where workflow Docker images are hosted."}
        zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
    }
}


task merge_and_plot_qc_metrics {
    input {
        String output_prefix
        Array[String] merge_keys
        Array[File] prepped_adata_objects
        String container_registry
        String zones
    }

    Int mem_gb = ceil(size(prepped_adata_objects, "GB") * 2.4 + 20)
    Int disk_size = ceil(size(prepped_adata_objects, "GB") * 3 + 50)

    command <<<
        set -euo pipefail

        paste ~{write_lines(merge_keys)} ~{write_lines(prepped_adata_objects)} > adata_samples_paths.tsv

        /usr/bin/time -v \
        merge_and_plot_qc \
            --adata-objects-fofn adata_samples_paths.tsv \
            --adata-output ~{output_prefix}.merged_unfiltered.h5ad \
            --output-metadata-file ~{output_prefix}.initial_metadata.csv

        mv "plots/violin_n_genes_by_counts.png" "plots/~{output_prefix}.n_genes_by_counts.violin.png"
        mv "plots/violin_total_counts.png" "plots/~{output_prefix}.total_counts.violin.png"
        mv "plots/violin_pct_counts_mt.png" "plots/~{output_prefix}.pct_counts_mt.violin.png"
        mv "plots/violin_pct_counts_rb.png" "plots/~{output_prefix}.pct_counts_rb.violin.png"
        mv "plots/violin_doublet_score.png" "plots/~{output_prefix}.doublet_score.violin.png"

    >>>

    output {
        File merged_adata_object = "~{output_prefix}.merged_unfiltered.h5ad"
        File qc_initial_metadata_csv = "~{output_prefix}.initial_metadata.csv"
        Array[File] qc_plots_png = [
            "plots/~{output_prefix}.n_genes_by_counts.violin.png",
            "plots/~{output_prefix}.total_counts.violin.png",
            "plots/~{output_prefix}.pct_counts_mt.violin.png",
            "plots/~{output_prefix}.pct_counts_rb.violin.png",
            "plots/~{output_prefix}.doublet_score.violin.png"
        ]
    }

    runtime {
        docker: "~{container_registry}/parsebio_tools:0.3.1"
        cpu: 4
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} SSD"
        bootDiskSizeGb: 40
        preemptible: 3
        zones: zones
    }

    meta {
        description: "Merges Parse Bio sample-level AnnData objects to a single AnnData object and generates pre-filtering QC violin plots for key metrics."
    }

    parameter_meta {
        output_prefix: {help: "Prefix for output files of the form <dataset_id>.<run_id>."}
        merge_keys: {help: "Concat key per AnnData object, in the same order as prepped_adata_objects (<run_id>.<parse_sample>); becomes the cell-ID suffix."}
        prepped_adata_objects: {help: "An array of Parse Sample-level AnnData objects to merge."}
        container_registry: {help: "Container registry where workflow Docker images are hosted."}
        zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
    }
}


task filter {
    input {
        String output_prefix
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

    Int calc_mem_gb = ceil(size(merged_adata_object, "GB") * 18 + 20)
    Int mem_gb = if calc_mem_gb > 624 then 624 else calc_mem_gb
    Int disk_size = ceil(size(merged_adata_object, "GB") * 4 + 20)

    command <<<
        set -euo pipefail

        /usr/bin/time -v \
        drop_unassigned_cells.py \
            --adata-input ~{merged_adata_object} \
            --adata-output ~{output_prefix}.assigned.h5ad \
            --output-removed-cells-file ~{output_prefix}.unassigned_removed_cells.csv

        /usr/bin/time -v \
        filter \
            --adata-input ~{output_prefix}.assigned.h5ad \
            --pct-counts-mt-max ~{pct_counts_mt_max} \
            --doublet-score-max ~{doublet_score_max} \
            --total-counts-limits ~{sep=' ' total_counts_limits} \
            --n-genes-by-counts-limits ~{sep=' ' n_genes_by_counts_limits} \
            --n-mads-lower ~{n_mads_lower} \
            --n-mads-upper ~{n_mads_upper} \
            --adata-output ~{output_prefix}.filtered.h5ad \
            --output-removed-cells-file ~{output_prefix}.filtered_out_cells.csv
    >>>

    output {
        File filtered_adata_object = "~{output_prefix}.filtered.h5ad"
        File unassigned_removed_cells_csv = "~{output_prefix}.unassigned_removed_cells.csv"
        File filtered_out_cells_csv = "~{output_prefix}.filtered_out_cells.csv"
    }

    runtime {
        docker: "~{container_registry}/parsebio_tools:0.3.1"
        cpu: 4
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
        output_prefix: {help: "Prefix for output files of the form <dataset_id>.<run_id>."}
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


task map_cell_types {
    input {
        String output_prefix
        File filtered_adata_object

        String mmc_taxonomy
        File allen_brain_mmc_precomputed_stats_h5
        File? allen_brain_mmc_marker_genes_json

        String container_registry
        String zones
    }

    String mmc_mapping_mode = if defined(allen_brain_mmc_marker_genes_json) then "mmc_markers_mapping" else "mmc_otf_mapping"
    String mmc_output_prefix = "~{output_prefix}.~{mmc_mapping_mode}.~{mmc_taxonomy}"

    # Marker lookup mode skips reference marker generation (the expensive bit). 
    # If not provided, on-the-fly marker generation will be used, which also
    # requires more scratch storage that scales with size of the reference taxonomy.
    Boolean use_marker_lookup = defined(allen_brain_mmc_marker_genes_json)
    
    Int threads = if use_marker_lookup then 4 else 16
    
    Int calc_mem_gb = if use_marker_lookup
        then ceil(size([filtered_adata_object, allen_brain_mmc_precomputed_stats_h5], "GB") * 3 + 16)
        else ceil(size([filtered_adata_object, allen_brain_mmc_precomputed_stats_h5], "GB") * 8 + 100)
    Int mem_gb = if calc_mem_gb > 624 then 624 else calc_mem_gb

    Int disk_size = if use_marker_lookup
        then ceil(size([filtered_adata_object, allen_brain_mmc_precomputed_stats_h5], "GB") * 3 + 20)
        else ceil(size([filtered_adata_object, allen_brain_mmc_precomputed_stats_h5], "GB") * 4 + 250)

    command <<<
        set -euo pipefail

        /usr/bin/time -v \
        mmc \
            --adata-input ~{filtered_adata_object} \
            --mmc-precomputed-stats ~{allen_brain_mmc_precomputed_stats_h5} \
            --n-processors ~{threads} \
            --output-prefix ~{mmc_output_prefix} \
            ~{if defined(allen_brain_mmc_marker_genes_json) then "--mmc-marker-genes " + allen_brain_mmc_marker_genes_json else ""}
    >>>

    output {
        File mmc_extended_results_json = "~{mmc_output_prefix}.extended_results.json"
        File mmc_results_csv = "~{mmc_output_prefix}.results.csv"
        File mmc_log_txt = "~{mmc_output_prefix}.log.txt"
    }

    runtime {
        docker: "~{container_registry}/parsebio_tools:0.3.1"
        cpu: threads
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 0  # TODO: update to 3 after cached
        bootDiskSizeGb: 40
        zones: zones
    }

    meta {
        description: "Assigns cell type labels with the Allen Brain Cell Atlas MapMyCells (MMC) tool, using on-the-fly markers by default or a precomputed marker lookup if provided. Writes the per-cell results CSV, the extended results JSON (including the marker genes used) and the log."
    }

    parameter_meta {
        output_prefix: {help: "Prefix for output files of the form <dataset_id>.<run_id>."}
        filtered_adata_object: {help: "QC-filtered AnnData object."}
        mmc_taxonomy: {help: "Cell type taxonomy of the precomputed stats reference; appended to MMC output filenames. Must match allen_brain_mmc_precomputed_stats_h5."}
        allen_brain_mmc_precomputed_stats_h5: {help: "A precomputed statistics file from the Allen Brain Cell Atlas containing reference statistics (the average gene expression profile per cell type cluster and cell type taxonomy)."}
        allen_brain_mmc_marker_genes_json: {help: "Optional marker gene lookup (parent node -> marker genes), e.g. extracted from a previous run's extended_results.json. If provided, skips on-the-fly marker computation."}
        container_registry: {help: "Container registry where workflow Docker images are hosted."}
        zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
    }
}


task add_mapped_cell_types {
    input {
        String output_prefix
        File filtered_adata_object
        File mmc_results_csv

        String container_registry
        String zones
    }

    # transcriptional_phenotype get its labelling logic from sc-rnaseq-wf's workflow name.
    # "pmdbs_sc_rnaseq" selects the human branch, currently the only supported in this workflow.
    String sc_tools_workflow_name = "pmdbs_sc_rnaseq"

    # TODO: revisit, seems generous
    Int calc_mem_gb = ceil(size(filtered_adata_object, "GB") * 18 + 20)
    Int mem_gb = if calc_mem_gb > 624 then 624 else calc_mem_gb
    Int disk_size = ceil(size(filtered_adata_object, "GB") * 4 + 20)

    command <<<
        set -euo pipefail

        /usr/bin/time -v \
        transcriptional_phenotype \
            --adata-input ~{filtered_adata_object} \
            --workflow-name ~{sc_tools_workflow_name} \
            --mmc-results ~{mmc_results_csv} \
            --adata-output ~{output_prefix}.mmc.h5ad \
            --output-cell-types-file ~{output_prefix}.mmc_results.parquet
    >>>

    output {
        File mmc_adata_object = "~{output_prefix}.mmc.h5ad"
        File mmc_results_parquet = "~{output_prefix}.mmc_results.parquet"
    }

    runtime {
        docker: "~{container_registry}/parsebio_tools:0.3.1"
        cpu: 4
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 3
        bootDiskSizeGb: 40
        zones: zones
    }

    meta {
        description: "Joins MMC cell type mapping results onto the filtered AnnData object and exports a Parquet summary of per-cell type assignments."
    }

    parameter_meta {
        output_prefix: {help: "Prefix for output files of the form <dataset_id>.<run_id>."}
        filtered_adata_object: {help: "QC-filtered AnnData object."}
        mmc_results_csv: {help: "Cell type mapping results CSV produced by map_cell_types."}
        container_registry: {help: "Container registry where workflow Docker images are hosted."}
        zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
    }
}