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

Users are suggested to use [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html) to install these dependencies. For example, the following shell command could install a workflow handling ONT R9 and R10 sequencing data in a **conda** environment named **snakemake-tutorial**.
```bash
conda env create --name snakemake-tutorial --file R9.environment.yaml
```
## Activate and exit the environment
To activate the environment 
  ```bash
  conda activate snakemake-tutorial
  ```
To exit the environment (after finishing the usage of the pipeline), just execute
  ```bash
  conda deactivate
  ```
## Run the pipeline
The whole pipeline contains two separate workflows to handle ONT R9 and R10 sequencing data respectively, i.e., **R9.snakemake.py** for **R9** data and **R10.snakemake.py** for **R10** data.
### Configure input parameters
```bash
Usage:
snakemake -s R9.snakemake.test.py --cores 8 --config Ref=refname Num=5 Vcf=path/to/VCF Refseq=path/to/refseq Outdir=path/to/outputdir Bam=path/to/sorted/bam Tombo_dir=path/to/tombo_processed/fast5 Subsample=2000
Ref=refname                             The value of #CHROM in vcf file, e.g., 'Ref=chr1'
Num=5                                   The number of bases up- and down-stream that are centered around the variation loci, default=5
Vcf=path/to/VCF                         The file path to vcf file, e.g., 'Vcf=/data/res/lofreq.vcf'
Refseq=path/to/refseq                   The file path to reference sequence, e.g., 'Refseq=/database/COVID-19.fa'
Outdir=path/to/outputdir                The file path storing the output results and intermediate files, e.g., 'Outdir=/data/res'
Bam=path/to/sorted/bam                  The file path to sorted bam files, e.g., 'Bam=/data/res/sorted.bam'
Tombo_dir=path/to/tombo_processed/fast5 The file path to tombo-resquiggled single fats5 files, e.g., 'Tombo_dir=/data/fast5'
Subsample=2000                          The number to subsample from reads covering variation loci, should be larger than 200, default=2000
```
The working directory contains files named `multi-DegePrime.yaml`, `multiPrime-original.yaml` and `multiPrime.yaml`. These are the central file in which all user settings, paramter values and path specifications are stored. `multi-DegePrime.yaml` employs DEGEPRIME-1.1.0 for maximum coverage degenerate primer design (MC-DPD), `multiPrime-orignal.yaml` and `multiPrime.yaml` use multiPrime-core.py for MC-DPD or MC-DPD with error. During a run, all steps of the pipeline will retrieve their paramter values from these file. It follows the yaml syntax (find more information about yaml and it's syntax [here](http://www.yaml.org/)) what makes it easy to read and edit. The main principles are:
  - everything that comes after a `#` symbol is considered as comment and will not be interpreted
  - paramters are given as key-value pair, with `key` being the name and `value` the value of any paramter

Before starting the pipeline, open the `multiPrime.yaml` configuration file and set all options according as required. This should at least include:
  - **name of the input directory** - where are your input fasta files stored
	-input_dir: ["abs_path_to_input_dir"]
  - **name of the output directory** - where should the pipeline store the output files (the direcotry is created if not existing)
	-results_dir: ["abs_path_to_results_dir"]
  - **name of the log directory** - where should the pipeline store the log files
	-log_dir: ["abs_path_to_log_dir"]
  - **name of the scripts directory** - where should the pipeline store the scripts files
	-scripts_dir: ["abs_path_to"]/multiPrime/scripts
  - **name(s) of your input samples** - please note: If your sample is named `sample1.fa` then `sample1` will be kept as naming scheme throughout the entire run to indicate output files that belong to this input file, e.g. the pipeline will create a file called `sample1.fa`. If you have multiple input files, just follow the given pattern with one sample name per line (and a dash that indicates another list item).
  - **identity** - threshold for classification. please note: If you set 1, multiPrime will design candidate primer pairs for each fasta in input files. Suggestion: 0.7-0.8. 
  - **others** - for more information on the parameters, please refer to the YAML file.

# Start a run

Once you set up your configuration file, running the pipeline locally on your computer is as easy as invoking:
  ```bash
  sh run.sh
  ```
  maximal coverage degenerate primer design (MC-DPD). The approach employed DegePrime to design degenerate primers for the target sequence.
  ```bash
  snakemake --configfile multi-DegePrime.yaml -s multi-DegePrime.py --cores 10 --resources disk_mb=80000
  ```
  maximal coverage degenerate primer design with errors tolerant (MC-EDPD) or MC-DPD. MultiPrime-orignal is capable of avoiding mismatches that occur at the 3'end position. The approach used in multiPrime-orignal.yaml depends on the value of the "variation" parameter.

  If "variation" is set to 0, then multiPrime uses the MC-DPD approach to design degenerate primers for the target sequence. In this approach, the primer sequences are designed with prefect match (0-mismatch).

  If "variation" is set to a value greater than 0, then multiPrime uses the MC-EDPD approach to design degenerate primers with errors (mismatches) tolerance (1-mismatch or 2-mismatches). In this approach, the primer sequences are allowed to contain a limited number of errors (mismatches), which increases the probability of finding suitable primer sequences for the target sequence.
  ```bash
  snakemake --configfile multiPrime-orignal.yaml -s multiPrime-orignal.py --cores 10 --resources disk_mb=80000
  ```
  multiPrime is similiar to multiPrime-orignal, but it enables easy avoidance of mismatches at any position.
  ```bash
  snakemake --configfile multiPrime.yaml -s multiPrime.py --cores 10 --resources disk_mb=80000
  ```

