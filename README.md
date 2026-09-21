# split-pipe-wf

Repo for testing and developing a workflow to process scRNA-seq data using Parse Bioscience's `split-pipe` pipeline.

Common workflows, tasks, utility scripts, and docker images reused across harmonized ASAP workflows are defined in [the wf-common repository](https://github.com/ASAP-CRN/wf-common).


# Table of contents

- [Workflows](#workflows)
- [Inputs](#inputs)
- [Outputs](#outputs)
    - [Output structure](#output-structure)
- [Docker images](#docker-images)


> **Status: work in progress.** The `mkref → all → comb` skeleton runs end to end,
> but outputs are not yet declared individually (see [Outputs](#outputs)), and no
> downstream QC, doublet detection, or bucket staging is implemented.


# Workflows

Workflows are defined in [the `workflows` directory](workflows).

Workflow name:
- `split_pipe_entrypoint`

**Entrypoint**: [workflows/main.wdl](workflows/main.wdl)

**Input template**: [workflows/inputs_template.json](workflows/inputs_template.json)

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


# Inputs

An input template file can be found at
[workflows/inputs_template.json](workflows/inputs_template.json).

| Type | Name | Description |
| :- | :- | :- |
| Array[[Sublibrary](#sublibrary)] | sublibraries | The set of sublibraries to process, and their associated reads. |
| String | organism | Organism; used to name the reference. Options: 'human' (hg38) or 'mouse' (mm10). |
| File? | genome_ref_tar_gz | Prebuilt `split-pipe` mkref index, as a gzipped tarball. If supplied, `make_reference` is skipped; if omitted, an index is built from `genome_fasta` and `genes_gtf`. |
| File | genome_fasta | Reference genome FASTA. Only used when `genome_ref_tar_gz` is omitted. |
| File | genes_gtf | Reference gene annotation GTF. Only used when `genome_ref_tar_gz` is omitted. |
| String | chemistry | Parse Biosciences chemistry version. Options: 'v1'–'v4'. |
| String | kit | Parse Biosciences kit. Options: 'WT_mini', 'WT', 'WT_mega', 'WT_mega_384', 'WT_penta', 'WT_penta_384'. |
| File | sample_list_tsv | Sample table; one `sample_name<whitespace>wells` line per sample, where wells are **round-1** wells: how samples were aliquoted onto the plate prior to the rounds of barcoding (e.g. `sample_1 A1-A12,B1-B3`). |
| File? | parameter_file | `split-pipe --parfile`; one `key value` pair per line. |
| String | container_registry | Container registry where workflow Docker images are hosted. |
| String? | zones | GCP zones where compute will take place. ['us-central1-c us-central1-f'] |

The sample table is a property of the kit run, not of any one sublibrary — every
sublibrary of a kit run takes the same file. 

## Structs

### Sublibrary

| Type | Name | Description |
| :- | :- | :- |
| String | sublibrary_id | Identifier for the sublibrary; used to name its output directory and tarball. |
| Array[File]+ | fastq_R1s | Read 1 (cDNA) FASTQs. Multiple files are accepted where a sublibrary was sequenced across lanes or runs. |
| Array[File]+ | fastq_R2s | Read 2 (barcode) FASTQs, in the same order as `fastq_R1s`. |


# Outputs

| Type | Name | Description |
| :- | :- | :- |
| File | combined_sublibs_tar_gz | The `split-pipe --mode comb` output directory, as a gzipped tarball. |

> **WIP:** the workflow currently emits a single tarball. Declaring the
> individual outputs — per-sample count matrices, the annotated BAMs, the
> AnnData objects and the HTML summaries — and staging them into raw and staging
> buckets is not yet implemented.


# Docker images

Docker images are defined in [the `docker` directory](docker). Each image must minimally define a `build.env` file and a `Dockerfile`.

Structure:
```bash
docker
└── split_pipe
    ├── build.env
    └── Dockerfile
```

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