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
The working directory contains files named `R9.snakemake.config.yaml` and `R10.snakemake.config.yaml`, which stores key input parameters for each workflow.
### Configure input parameters for R9 workflow
There are two ways to configure input parameters for R9 workflow.

(1) Via shell command line

Users could define customized input paramaters using **--config** option in Snakemake command line.
```bash
Usage:
snakemake -s R9.snakemake.py --cores 8 --config Ref=refname Num=5 Vcf=path/to/VCF Refseq=path/to/refseq Outdir=path/to/outputdir Bam=path/to/sorted/bam Tombo_dir=path/to/tombo_processed/fast5 Subsample=2000
Ref=refname                             The value of #CHROM in vcf file, e.g., 'Ref=chr1'
Num=5                                   The number of bases up- and down-stream that are centered around the variation loci, default=5
Vcf=path/to/VCF                         The file path to vcf file, e.g., 'Vcf=/data/res/lofreq.vcf'
Refseq=path/to/refseq                   The file path to reference sequence, e.g., 'Refseq=/database/COVID-19.fa'
Outdir=path/to/outputdir                The file path storing the output results and intermediate files, e.g., 'Outdir=/data/res'
Bam=path/to/sorted/bam                  The file path to sorted bam files, e.g., 'Bam=/data/res/sorted.bam'
Tombo_dir=path/to/tombo_processed/fast5 The file path to tombo-resquiggled single fats5 files, e.g., 'Tombo_dir=/data/fast5'
Subsample=2000                          The number to subsample from reads covering variation loci, should be larger than 200, default=2000
```
(2) Edit config.yaml

Users could also define customized input paramaters by editing R9.snakemake.config.yaml.
```yaml
Bam: "/public/data1/yefq/data/fast5/20220703_WGA_twist/processed/20230426_Guppy621_comparison/Sce20_guppy_sup_aligned.softclip_trimmed.endtrim10_minimap2_align.mapped.sorted.bam"
Ref: "Zymo_Saccharomyces_cerevisiae_Seq5_ref"
Num: "5"
Tombo_dir: "/public/data1/yefq/data/fast5/20220703_WGA_twist/processed/20230426_guppy_sup_basecalled/Sce20/workspace/fast5_pass_single/all_single_fast5s"
Subsample: "2000"
Vcf: "/public/data1/yefq/data/fast5/20220703_WGA_twist/processed/20230426_Guppy621_comparison/Sce20_guppy_sup_aligned.test.vcf" #20240226 updated
Refseq: "/public/data1/yefq/data/Refs/Zymo_Saccharomyces_cerevisiae_Seq5_ref.fa" #240225 updated
Outdir: "/public/data1/yefq/data/fast5/20220703_WGA_twist/processed/20230426_Guppy621_comparison/snakemake-tutorial/data/test" #240225 updated
```
Users should note that, **config values can be overwritten via the command line** even when it has deen defined in the config.yaml.
### Configure input parameters for R10 workflow
There are also two ways to configure input parameters for R10 workflow.

(1) Via shell command line

Users could define customized input paramaters using **--config** option in Snakemake command line.
```bash
Usage:
snakemake -s R10.snakemake.py --cores n --config Ref=refname Num=5 Vcf=path/to/VCF Refseq=path/to/refseq Outdir=path/to/outputdir Bam=path/to/sorted/bam
Ref=refname                   The value of #CHROM in vcf file, e.g., 'Ref=chr1'
Num=5                         The number of bases up- and down-stream that are centered around the variation loci, default=5
Vcf=path/to/VCF               The file path to vcf file, e.g., 'Vcf=/data/res/lofreq.vcf'
Refseq=path/to/refseq         The file path to reference sequence, e.g., 'Refseq=/database/COVID-19.fa'
Outdir=path/to/outputdir      The file path storing the output results and intermediate files, e.g., 'Outdir=/data/res'
Bam=path/to/sorted/bam        The file path to sorted bam files, e.g., 'Bam=/data/res/sorted.bam'
```
(2) Edit config.yaml

Users could also define customized input paramaters by editing R10.snakemake.config.yaml.
```yaml
Ref: "Zymo_Saccharomyces_cerevisiae_Seq5_ref"
Num: "5"
Vcf: "/public/data1/yefq/data/fast5/20220703_WGA_twist/processed/20230426_Guppy621_comparison/Sce20_guppy_sup_aligned.test.vcf" #20240226 updated
Refseq: "/public/data1/yefq/data/Refs/Zymo_Saccharomyces_cerevisiae_Seq5_ref.fa" #240225 updated
Outdir: "/public/data1/yefq/data/fast5/20220703_WGA_twist/processed/20230426_Guppy621_comparison/snakemake-tutorial-R10/data/test" #240225 updated
Bam: "/public/data1/yefq/data/fast5/20220703_WGA_twist/processed/20230426_Guppy621_comparison/Sce20_guppy_sup_aligned.softclip_trimmed.endtrim10_minimap2_align.mapped.sorted.bam"
```
Users should note that, **config values can be overwritten via the command line** even when it has deen defined in the config.yaml.

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

