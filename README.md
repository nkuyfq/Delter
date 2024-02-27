# ONT_Deletion_Filter
A tool to filter short artificial deletion variations by Oxford Nanopore Technologies (ONT) R9 and R10 flow cells and chemistries.
## Requirements
The tool has been tested on Ubuntu 20.04 with 256GB RAM, 64 CPU cores and a NVIDIA GPU with 48GB RAM. The minimal requirements should be >= 64GB RAM and a NVIDIA GPU with >= 8GB RAM. Other operating systems like Windows or Mac were not tested.
ONT softwares like Guppy, Tombo, ont-fast5-api should be pre-installed before generating Tombo-resquiggled single-read fast5 files.
## Installation
The tool runs via Snakemake workflows. Users must install workflow dependencies including [Snakemake](https://snakemake.readthedocs.io/en/latest/tutorial/tutorial.html) before using the pipeline. The workflow dependencies, which stored in a file named environment.yaml, are listed as below:

channels:
  - conda-forge
  - bioconda
  - anaconda

dependencies:
  - snakemake-minimal >=7.3
  - graphviz
  - seaborn
  - numpy
  - pandas
  - h5py
  - scipy
  - samtools =1.15
  - r-essentials
  - r-base
  - bioconductor-shortread
  - r-stringr
  - r-dplyr
  - r-vegan 

Users are suggested to use [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html) to install these dependencies. For example, the following shell command could install a workflow handling ONT R9 and R10 sequencing data.
```bash
conda env create --name snakemake-tutorial --file R9.environment.yaml
```

## Installation
