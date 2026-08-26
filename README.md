# split-pipe-wf

Repo for testing and developing a workflow to process scRNA-seq data using Parse Bioscience's `split-pipe` pipeline.

Common workflows, tasks, utility scripts, and docker images reused across harmonized ASAP workflows are defined in [the wf-common repository](https://github.com/ASAP-CRN/wf-common).


# Table of contents

- [Workflows](#workflows)
- [Inputs](#inputs)
- [Outputs](#outputs)
    - [Output structure](#output-structure)
- [Docker images](#docker-images)


# Workflows

TODO

# Inputs

TODO

# Outputs

TODO

# Docker images

Docker images are defined in [the `docker` directory](docker). Each image must minimally define a `build.env` file and a `Dockerfile`.

Structure:
```bash
docker
└── split_pipe
    ├── build.env
    ├── Dockerfile
    └── fetch_split_pipe.sh

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
./build_docker_images -d docker/split-pipe

# Build all images in the `docker` directory
./build_docker_images -d docker

# Build and push all images in the docker directory, using the `dnastack` container registry
./build_docker_images -d docker -c dnastack -p
```


## Tool and library versions

| Image | Major tool versions | Links |
| :- | :- | :- |
| split-pipe | <ul><li>[split-pipe v1.8.2 (login required)](https://support.parsebiosciences.com/hc/en-us/articles/17200056667924-Pipeline-Download-Current-Version)</li><li>[google-cloud-cli 572.0.0](https://docs.cloud.google.com/sdk/docs/release-notes#57200_2026-06-09)</li><li>[python 3.12.8](https://www.python.org/downloads/release/python-3128/)</li><li>[miniconda 26.5.3](https://docs.anaconda.com/miniconda/miniconda-release-notes/)</li></ul> | [Dockerfile](https://github.com/ASAP-CRN/split-pipe-wf/tree/main/docker/split_pipe) |