# split-pipe-wf

Repo for testing and developing a workflow to process scRNA-seq data using Parse Bioscience's `split-pipe` pipeline.

Common workflows, tasks, utility scripts, and docker images reused across harmonized ASAP workflows are defined in [the wf-common repository](https://github.com/ASAP-CRN/wf-common).


# Table of contents

- [Workflows](#workflows)
- [Inputs](#inputs)
    - [Postprocess inputs](#postprocess-inputs)
    - [Structs](#structs)
- [Outputs](#outputs)
    - [Postprocess outputs](#postprocess-outputs)
- [Docker images](#docker-images)


> **Status: work in progress.** `main.wdl` runs every Parse run from FASTQs to a filtered,
> cell-typed AnnData object and merges the runs into one object. Declaring the remaining
> outputs (count matrices, BAMs, reports) and staging them into raw and staging buckets is
> not yet implemented.


# Workflows

Workflows are defined in [the `workflows` directory](workflows).

Workflow names:
- `split_pipe_entrypoint` ([workflows/main.wdl](workflows/main.wdl)): the full workflow. Builds
  the reference once, then for each Parse run calls `split_pipe` and `postprocess`, and finally
  merges the runs (`merge_runs`).
- `split_pipe` ([workflows/split_pipe/split_pipe.wdl](workflows/split_pipe/split_pipe.wdl)):
  `split-pipe` processing for one Parse run, from FASTQs to the combined (`comb`) output.
- `postprocess` ([workflows/postprocess/postprocess.wdl](workflows/postprocess/postprocess.wdl)):
  from one run's `comb` output to a filtered, cell-typed AnnData object.

`split_pipe` and `postprocess` can also be run on their own.

**Entrypoint**: [workflows/main.wdl](workflows/main.wdl)

**Input template**: [workflows/inputs_template.json](workflows/inputs_template.json)

The Python scripts used by postprocessing are in [docker/parsebio_tools/scripts](docker/parsebio_tools/scripts). Merging, QC plotting and filtering reuse the `sc_tools` scripts from [sc-rnaseq-wf](https://github.com/ASAP-CRN/sc-rnaseq-wf), so QC metric definitions and filter thresholds match the 10X workflow.

## split-pipe processing

A Parse experiment is sequenced as several **sublibraries**, each an aliquot of
one mixed pool of barcoded nuclei. Each is processed independently and the
results are combined afterwards:

```
make_reference           once per reference; builds the STAR index (mkref)
        │
        ├── split_pipe   scattered, one call per sublibrary (mode: all)
        │
combine_sublibraries     gather; resolves cells to samples across sublibraries (mode: comb)
```

`make_reference` is skipped when a prebuilt index is supplied via
`genome_ref_tar_gz`.

`combine_sublibraries` requires at least two sublibraries, and any subset of a kit run may be
combined: each sublibrary contains a proportion of cells from all samples.

## Postprocessing

Run once per Parse run (one kit's worth of samples, e.g. `PB-001`), on that run's `comb` output.

```
prep_parse_sample          scattered, one call per Parse sample (e.g. pool_1)
  ├── assign_donors.py         cells → subject, from a donor demultiplexing table (vireo)
  │                            or a single known subject (un-pooled sample)
  └── prep_parse_sample.py     AnnData from DGE_filtered; QC metrics, scrublet scores,
                               ASAP SAMPLE/SUBJECT metadata
        │
merge_and_plot_qc_metrics  gather; concatenates Parse samples, pre-filter QC violin plots
        │                  (sc_tools merge_and_plot_qc)
filter
  ├── drop_unassigned_cells.py drops cells without an ASAP sample (vireo doublet/unassigned)
  └── filter                   sc_tools filter: mito, doublet score, count/gene limits and
                               per-sample MAD outliers
        │
map_cell_types             MapMyCells (cell_type_mapper), on-the-fly markers from gene symbols
        │                  (sc_tools mmc)
add_mapped_cell_types      attaches MMC labels and phenotypes (sc_tools transcriptional_phenotype)
```

After every run is postprocessed, `merge_runs` (in `main.wdl`) concatenates the per-run
cell-typed objects into one object (`merge_runs.py`). It checks that runs share the gene set,
carry distinct `batch` values and have unique cell IDs. Runs can be merged at the end because
no sample spans two runs and MMC labels each cell independently.

Notes:

- **Parse sample ≠ ASAP sample.** A Parse sample is what was loaded into round-1 wells. In
  pooled designs it holds several donors, which are resolved per cell by `assign_donors.py`.
  The ASAP sample ID is then looked up from (batch, subject); `split-pipe`'s own `sample`
  column is kept as `parse_sample`.
- **Cells come from `DGE_filtered`**, i.e. `split-pipe`'s own cell calling, for every Parse
  sample, pooled or not.
- **Doublets.** Scrublet runs on every called cell, including vireo doublets, so its
  `doublet_score` sees the full doublet population. Vireo detects doublets between donors,
  and scrublet detects doublets between cell types; both are applied.
- **Cell IDs** are `<bc_wells>_<run_id>.<parse_sample>_<n>` (e.g. `01_01_15__s1_PB-001.pool_1_0`),
  unique across Parse runs. The trailing `_<n>` is the merge chunk added by sc_tools
  `merge_and_plot_qc`, as in sc-rnaseq-wf.


# Inputs

An input template file can be found at
[workflows/inputs_template.json](workflows/inputs_template.json).

The reference is built from the Cell Ranger GRCh38-2024-A FASTA and GTF
(`fasta/genome.fa`, `genes/genes.gtf.gz`), the same annotation sc-rnaseq-wf uses, so both
workflows count the same gene set.

Inputs that differ between Parse runs (kit, chemistry, batch, sublibraries, Parse samples)
are fields of each [ParseRun](#parserun); everything else is shared.

| Type | Name | Description |
| :- | :- | :- |
| String | dataset_id | The unique name used to identify the dataset across the CRN Cloud. |
| String | team_id | ASAP team identifier; used in `batch_id`. |
| Array[[ParseRun](#parserun)] | parse_runs | The Parse runs (kits) in the dataset. |
| String | genome_name | Genome name; used to name the reference. Options: 'hg38' (human) or 'mm10' (mouse). |
| File? | genome_ref_tar_gz | Prebuilt `split-pipe` mkref index, as a gzipped tarball. If supplied, `make_reference` is skipped; if omitted, an index is built once from `genome_fasta` and `genes_gtf` and shared by all runs. |
| File | genome_fasta | Reference genome FASTA. Only used when `genome_ref_tar_gz` is omitted. |
| File | genes_gtf | Reference gene annotation GTF. Only used when `genome_ref_tar_gz` is omitted. |
| File? | parameter_file | `split-pipe --parfile`; one `key value` pair per line. |
| File | sample_csv | ASAP `SAMPLE.csv`. |
| File | subject_csv | ASAP `SUBJECT.csv`. |
| Int | pct_counts_mt_max | Maximum percentage of mitochondrial gene counts allowed per cell. [10] |
| Float | doublet_score_max | Maximum scrublet doublet score allowed per cell. [0.2] |
| Array[Int] | total_counts_limits | Absolute minimum and maximum total UMI counts per cell. [500, 150000] |
| Array[Int] | n_genes_by_counts_limits | Absolute minimum and maximum number of genes detected per cell. [300, 15000] |
| Float | n_mads_lower | Number of MADs below the per-sample median allowed for total counts and genes detected. [3] |
| Float | n_mads_upper | Number of MADs above the per-sample median allowed for total counts and genes detected. [5] |
| String | mmc_taxonomy | Taxonomy of the MapMyCells reference; appended to MMC output names (e.g. `siletti_20240330`). |
| File | allen_brain_mmc_precomputed_stats_h5 | MapMyCells precomputed statistics for the reference taxonomy. |
| File? | allen_brain_mmc_marker_genes_json | Marker gene lookup; if supplied, skips on-the-fly marker computation. |
| String | container_registry | Container registry where workflow Docker images are hosted. |
| String? | zones | GCP zones where compute will take place. ['us-central1-c us-central1-f'] |

## Postprocess inputs

| Type | Name | Description |
| :- | :- | :- |
| File | combined_sublibs_tar_gz | The `split-pipe --mode comb` output directory for one Parse run, as a gzipped tarball. |
| String | dataset_id | ASAP dataset identifier; used in output names and `batch_id`. |
| String | team_id | ASAP team identifier; used in `batch_id`. |
| String | run_id | Parse run ID (e.g. `PB-001`); used in output names and cell IDs. |
| String | batch | `SAMPLE.csv` batch value for this Parse run (e.g. `Batch1`). |
| Array[[ParseSample](#parsesample)] | parse_samples | The Parse samples to process, and how to assign their cells to donors. |
| File | sample_csv | ASAP `SAMPLE.csv`. |
| File | subject_csv | ASAP `SUBJECT.csv`. |
| File? | sample_map | CSV of `parse_sample,subject_id,sample_id`. Required when a subject has more than one sample in the batch, e.g. two brain regions loaded in different Parse samples. See [Sample map](#sample-map). |
| Int | pct_counts_mt_max | Maximum percentage of mitochondrial gene counts allowed per cell. [10] |
| Float | doublet_score_max | Maximum scrublet doublet score allowed per cell. [0.2] |
| Array[Int] | total_counts_limits | Absolute minimum and maximum total UMI counts per cell, applied on top of the MAD-based thresholds. [500, 150000] |
| Array[Int] | n_genes_by_counts_limits | Absolute minimum and maximum number of genes detected per cell, applied on top of the MAD-based thresholds. [300, 15000] |
| Float | n_mads_lower | Number of MADs below the per-sample median allowed for total counts and genes detected. [3] |
| Float | n_mads_upper | Number of MADs above the per-sample median allowed for total counts and genes detected. [5] |
| String | mmc_taxonomy | Taxonomy of the MapMyCells reference; appended to MMC output names. |
| File | allen_brain_mmc_precomputed_stats_h5 | MapMyCells precomputed statistics for the reference taxonomy. |
| File? | allen_brain_mmc_marker_genes_json | Marker gene lookup; if supplied, skips on-the-fly marker computation. |
| String | container_registry | Container registry where workflow Docker images are hosted. |
| String? | zones | GCP zones where compute will take place. ['us-central1-c us-central1-f'] |

Filter and MMC defaults match sc-rnaseq-wf.

### Sample map

Vireo assigns each cell to a *subject* by genotype, but ASAP metadata is keyed by *sample*. When a subject contributed more than one sample to the same Parse run (e.g. substantia nigra and cingulate cortex), genotype alone cannot tell the two samples apart. What does tell them apart is which Parse sample (set of round-1 wells) the cell was loaded in. `sample_map` is a curated input that records that link: one row per ASAP sample as (Parse sample, ASAP subject, ASAP sample). Cells are then joined to `SAMPLE.csv` on (batch, Parse sample, subject). `SAMPLE.csv` cannot provide the link itself, because its `pool_id` is the sublibrary and every sublibrary contains cells from every sample. Without a map, prep fails and names the ambiguous subjects.

## Structs

### ParseRun

One Parse run: one kit's worth of samples, processed by `split_pipe` and `postprocess`.

| Type | Name | Description |
| :- | :- | :- |
| String | run_id | Parse run ID (e.g. `PB-001`); used in output names and cell IDs. |
| String | kit | Parse Biosciences kit. Options: 'WT_mini', 'WT', 'WT_mega', 'WT_mega_384', 'WT_penta', 'WT_penta_384'. |
| String | chemistry | Parse Biosciences chemistry version. Options: 'v1'–'v4'. |
| String | batch | `SAMPLE.csv` batch value for this run (e.g. `Batch1`). |
| File | sample_list_tsv | Sample table; one `sample_name<whitespace>wells` line per Parse sample, where wells are **round-1** wells: how samples were aliquoted onto the plate prior to the rounds of barcoding (e.g. `sample_1 A1-A12,B1-B3`). Every sublibrary of the run takes the same file. |
| File? | sample_map | See [Sample map](#sample-map). |
| Array[[Sublibrary](#sublibrary)] | sublibraries | The run's sublibraries and their reads. |
| Array[[ParseSample](#parsesample)] | parse_samples | The run's Parse samples, and how to assign their cells to donors. |

### Sublibrary

| Type | Name | Description |
| :- | :- | :- |
| String | sublibrary_id | Identifier for the sublibrary; used to name its output directory and tarball. |
| Array[File]+ | fastq_R1s | Read 1 (cDNA) FASTQs. Multiple files are accepted where a sublibrary was sequenced across lanes or runs. |
| Array[File]+ | fastq_R2s | Read 2 (barcode) FASTQs, in the same order as `fastq_R1s`. |

### ParseSample

Exactly one of `vireo_donor_ids` and `subject_id` must be set.

| Type | Name | Description |
| :- | :- | :- |
| String | parse_sample_name | Parse sample name, as in the sample table; must match a `comb` output subdirectory. |
| File? | vireo_donor_ids | Vireo donor assignment TSV for a pooled Parse sample (columns `cell`, `donor_id`, ...); `cell` must match `split-pipe`'s `bc_wells`, including the `__s<n>` sublibrary suffix. |
| String? | subject_id | Subject ID for an un-pooled Parse sample; assigned to every cell. |

`team-voet-pmdbs-sn-rnaseq-parsebio` vireo notes:
> `__s<n>` suffixes are assigned by `comb` from the order of the supplied sublibraries. The 
> provided vireo files to link cells to donors depend on this ordering: keep
> sublibraries in the same order as the run the vireo table was computed on;
> `assign_donors.py` fails if fewer than 10% of called cells match the vireo table (a wrong
> order matches ~1–2%), and warns below 50%. A partial match is expected where `split-pipe`
> calls more cells than vireo was run on.
>
> Vireo names a donor it learned without a supplied genotype by position (`donor<N>`).
> These are renamed to a subject only where the contributors confirmed the identity
> (`VIREO_DONOR_RENAMES` in `assign_donors.py`, keyed by dataset and Parse sample; currently
> Team Voet PB-002 `pool_5`, `donor9` → `ASA_137`). Any other `donor<N>` fails
> `assign_donors.py`.


# Outputs

| Type | Name | Description |
| :- | :- | :- |
| File | final_merged_adata | `${dataset_id}.merged.mmc.h5ad`; all runs' filtered, cell-typed cells, raw counts, all genes (see [AnnData contents](#anndata-contents)). |

> **WIP:** the workflow currently emits only the merged object; per-run outputs remain in the
> execution directory. Declaring the individual outputs — per-sample count matrices, the
> annotated BAMs, the HTML summaries and the QC reports — and staging them into raw and
> staging buckets is not yet implemented.

## Postprocess outputs

| Type | Name | Description |
| :- | :- | :- |
| Array[File] | prepped_adata | Per Parse sample: `${dataset_id}.${run_id}.${parse_sample}.prepped.h5ad`, all called cells with donor assignments, QC metrics, scrublet scores and ASAP metadata. |
| File | merged_adata | `${dataset_id}.${run_id}.merged_unfiltered.h5ad`; all Parse samples concatenated, before any filtering. |
| File | initial_metadata | `${dataset_id}.${run_id}.initial_metadata.csv`; `obs` of the merged object. |
| Array[File] | qc_plots | Pre-filter violin plots of `n_genes_by_counts`, `total_counts`, `pct_counts_mt`, `pct_counts_rb` and `doublet_score`. |
| File | filtered_adata_object | `${dataset_id}.${run_id}.filtered.h5ad`; cells passing donor assignment and all QC filters. |
| File | unassigned_removed_cells_csv | Cells removed for having no ASAP sample (vireo `doublet`/`unassigned`). |
| File | filtered_out_cells_csv | Cells removed by `filter`, with a column per filter they failed. Same layout as `unassigned_removed_cells_csv`. |
| File | mmc_results_csv, mmc_extended_results_json, mmc_log_txt | MapMyCells outputs: per-cell assignments, full results (incl. the marker genes used) and log. |
| File | mmc_adata_object | `${dataset_id}.${run_id}.mmc.h5ad`; the filtered object with cell type labels. |
| File | mmc_results_parquet | Per-cell cell type labels. |

### AnnData contents

Raw counts are in `X`. All genes are kept; no gene filtering is applied.

| Location | Contents |
| :- | :- |
| `var` index | Gene symbol, made unique with a `-<n>` suffix as in sc-rnaseq-wf |
| `var["gene_id"]`, `var["gene_id_versioned"]` | Ensembl ID without and with the GENCODE version suffix |
| `obs["sample"]` | ASAP sample ID (`ASAP_sample_id`); `NA` for vireo doublet/unassigned cells |
| `obs["parse_sample"]` | Parse sample the cell was loaded in (e.g. `pool_1`) |
| `obs["subject_id"]` | Assigned subject, or `doublet`/`unassigned` |
| `obs["vireo_*"]` | Vireo statistics (`prob_max`, `prob_doublet`, `n_vars` = informative SNPs, ...); `NA` for un-pooled samples |
| `obs["n_genes_by_counts"]`, `obs["total_counts"]`, `obs["pct_counts_mt"]`, `obs["pct_counts_rb"]` | QC metrics, as in sc-rnaseq-wf |
| `obs["doublet_score"]` | Scrublet doublet score, per Parse sample |
| `obs["mread_count"]` | Sequencing reads per cell (`split-pipe`) |
| `obs["bc1_well"]`, `obs["bc2_well"]`, `obs["bc3_well"]` | Plate wells for each barcoding round |
| `obs["sex"]`, `obs["brain_region_level_1"]`–`obs["brain_region_level_3"]` | From ASAP metadata, lower-cased |
| `obs["batch"]`, `obs["team"]`, `obs["dataset"]`, `obs["batch_id"]` | Run-level metadata; `batch_id` = `${team}_${dataset}_${batch}`, as in sc-rnaseq-wf |
| `obs["supercluster_name"]`, `obs["cluster_name"]`, `obs["subcluster_name"]`, `obs["rho"]`, `obs["prob"]` | MapMyCells assignment and confidence (cell-typed objects only) |
| `obs["cell_type"]`, `obs["phenotype"]` | Labels from sc_tools `transcriptional_phenotype`, as in sc-rnaseq-wf (cell-typed objects only) |


# Docker images

Docker images are defined in [the `docker` directory](docker). Each image must minimally define a `build.env` file and a `Dockerfile`.

Structure:
```bash
docker
├── parsebio_tools
│   ├── build.env
│   ├── Dockerfile
│   └── scripts
└── split_pipe
    ├── build.env
    └── Dockerfile
```

`parsebio_tools` extends the sc-rnaseq-wf `sc_tools` image with the postprocessing scripts. The base is pinned by digest in `build.env` (`SC_TOOLS_DIGEST`), so library versions match the 10X workflow exactly.

## Acquiring split-pipe software

`split-pipe` is a proprietary pipeline owned by Parse Biosciences and licensed to DNAStack. Community users therefore must acquire their own license directly from Parse Biosciences if they wish to replicate this workflow. For DNAstack internal use, `split-pipe` is stored in a private Google bucket and must be saved into `/docker/split_pipe/split-pipe-{version_number}.zip`.


## The `build.env` file

Each target image is defined using the `build.env` file, which is used to specify the name and version tag for the corresponding Docker image. It must contain at minimum the following variables:

- `IMAGE_NAME`
- `IMAGE_TAG`

All variables defined in the `build.env` file will be made available as build arguments during Docker image build.

The `DOCKERFILE` variable may be used to specify the path to a Dockerfile if that file is not found alongside the `build.env` file, for example when multiple images use the same base Dockerfile definition.

## Building Docker images

Docker images can be build using the [`build_docker_images`](https://github.com/DNAstack/bioinformatics-scripts/blob/main/scripts/build_docker_images) script.

```bash
# Build a single image
./build_docker_images -d docker/split_pipe

# Build all images in the `docker` directory
./build_docker_images -d docker

# Build and push all images in the docker directory, using the `dnastack` container registry
./build_docker_images -d docker -c dnastack -p
```


## Tool and library versions

| Image | Major tool versions | Links |
| :- | :- | :- |
| split_pipe | <ul><li>[split-pipe v1.8.2 (login required)](https://support.parsebiosciences.com/hc/en-us/articles/17200056667924-Pipeline-Download-Current-Version)</li><li>[google-cloud-cli 572.0.0](https://docs.cloud.google.com/sdk/docs/release-notes#57200_2026-06-09)</li><li>[python 3.12.8](https://www.python.org/downloads/release/python-3128/)</li><li>[miniconda 26.5.3](https://docs.anaconda.com/miniconda/miniconda-release-notes/)</li></ul> | [Dockerfile](https://github.com/ASAP-CRN/split-pipe-wf/tree/main/docker/split_pipe) |
| parsebio_tools | <ul><li>Base: `sc_tools` 1.3.0, pinned by digest (scanpy, anndata and scrublet versions inherited)</li></ul> | [Dockerfile](https://github.com/ASAP-CRN/split-pipe-wf/tree/main/docker/parsebio_tools) |
