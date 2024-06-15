# ONT_Deletion_Filter-Delter
A tool to filter short artificial deletion variations by Oxford Nanopore Technologies (ONT) R9 and R10 flow cells and chemistries.
## Requirements
The tool has been tested on Ubuntu 20.04 with 256GB RAM, 64 CPU cores and a NVIDIA GPU with 48GB RAM. The minimal requirements should be >= 64GB RAM and a NVIDIA GPU with >= 8GB RAM. Other operating systems like Windows or Mac were not tested.

ONT softwares like [Guppy](https://community.nanoporetech.com/downloads), [Tombo](https://github.com/nanoporetech/tombo), and [ont-fast5-api](https://github.com/nanoporetech/ont_fast5_api) should be pre-installed before generating Tombo-resquiggled single-read fast5 files.
Users might run following commands to preprocess R9 fast5 files in shell terminal before running our pipeline. As these steps below need GPU support and might take a long time, our pipeline doesn't contain them.
```bash
#===basecalling the fast5 files===
ont-guppy/bin/guppy_basecaller -c ont-guppy/data/dna_r9.4.1_450bps_sup.cfg -i $fast5dir/barcode${barcode} -s guppy_sup_basecalled/barcode${barcode} -r --compress_fastq -x cuda:1,2 --gpu_runners_per_device 4 --chunks_per_runner 256 --num_callers 3 --fast5_out

#===preprocessing R9 fast5 files===
c=$(ls *.fast5 | wc -l)
declare -i count=$c-1
#multiread fast5 to single read fast5
multi_to_single_fast5 -i guppy_sup_basecalled/barcode${barcode}/workspace -s fast5_pass_single --threads 24 
#copy to new directory
cd fast5_pass_single
mkdir all_single_fast5s
for ((j=0;j<=$count;j=j+1))
do
  echo $j
  cp -r ./$j/*.fast5 all_single_fast5s
  rm -rf ./$j
done
#align to reference genome via tombo resquiggle
tombo resquiggle guppy_sup_basecalled/barcode${barcode}/workspace/fast5_pass_single/all_single_fast5s data/Refs/$refseq --processes 24 --overwrite --num-most-common-errors 5  --failed-reads-filename tombo_resquiggle_failed_fast5.txt
```
The fast5 files in the directory named **all_single_fast5s** could be employed in downstream workflow, which equals the input parameter **Tombo_dir** in shell command line or config yaml (details listed in **Configure input parameters** section). 
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

Users are suggested to use [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html) to install these dependencies. For example, the following shell command could install a workflow handling ONT R9 and R10 sequencing data in a **conda** environment named **Delter**.
```bash
conda env create --name Delter --file environment.yaml
```
## Activate and exit the environment
To activate the environment 
  ```bash
  conda activate Delter
  ```
To exit the environment (after finishing the usage of the pipeline), just execute
  ```bash
  conda deactivate
  ```
## Run the pipeline
The whole pipeline could handle ONT R9 and R10 sequencing data. The working directory contains file named `Delter.config.yaml`, which stores key input parameters for the workflow.
### Configure input parameters for the workflow
There are two ways to configure input parameters for this workflow.

(1) Via shell command line

Users could define customized input paramaters using **--config** option in Snakemake command line.
```bash
Usage:
snakemake -s Delter.py --cores 8 --config Ref=refname Num=5 Vcf=path/to/VCF Refseq=path/to/refseq Outdir=path/to/outputdir Bam=path/to/sorted/bam Tombo_dir=path/to/tombo_processed/fast5 Subsample=2000 Flowcell=R9 Strategy=Direct
Ref=refname                             The value of #CHROM in vcf file, e.g., 'Ref=chr1'
Num=5                                   The number of bases up- and down-stream that are centered around the variation loci, default=5
Vcf=path/to/VCF                         The file path to vcf file, e.g., 'Vcf=/data/res/lofreq.vcf'
Refseq=path/to/refseq                   The file path to reference sequence, e.g., 'Refseq=/database/COVID-19.fa'
Outdir=path/to/outputdir                The file path storing the output results and intermediate files, e.g., 'Outdir=/data/res'
Bam=path/to/sorted/bam                  The file path to sorted bam files, e.g., 'Bam=/data/res/sorted.bam'
Tombo_dir=path/to/tombo_processed/fast5 The file path to tombo-resquiggled single fats5 files, e.g., 'Tombo_dir=/data/fast5'
Subsample=2000                          The number to subsample from reads covering variation loci, should be larger than 200, default=2000
Flowcell=R9                             The version of flow cell, should be R9 or R10, default=R9
Strategy=Amplicon                       The sequencing strategy, should be Amplicon or Direct, default=Amplicon
```
(2) Edit config.yaml

Users could also define customized input paramaters by editing config.yaml.
```yaml
Bam: "/public/data1/yefq/data/fast5/20220703_WGA_twist/processed/20230426_Guppy621_comparison/Sce20_guppy_sup_aligned.softclip_trimmed.endtrim10_minimap2_align.mapped.sorted.bam"
Ref: "Zymo_Saccharomyces_cerevisiae_Seq5_ref"
Num: "5"
Tombo_dir: "/public/data1/yefq/data/fast5/20220703_WGA_twist/processed/20230426_guppy_sup_basecalled/Sce20/workspace/fast5_pass_single/all_single_fast5s"
Subsample: "2000"
Vcf: "/public/data1/yefq/data/fast5/20220703_WGA_twist/processed/20230426_Guppy621_comparison/Sce20_guppy_sup_aligned.test.vcf" 
Refseq: "/public/data1/yefq/data/Refs/Zymo_Saccharomyces_cerevisiae_Seq5_ref.fa" 
Outdir: "/public/data1/yefq/data/fast5/20220703_WGA_twist/processed/20230426_Guppy621_comparison/snakemake-tutorial/data/test" 
Flowcell: "R9"
Strategy: "Direct"
```
Users should note that, **config values can be overwritten via the command line** even when it has deen defined in the config.yaml.
### Start a run
Once the work directory and configuration files are set up, users can run the pipeline as easy as invoking:
```bash
cd work_dir
conda activate Delter
snakemake -s Delter.py --cores 8
```
Other Snakemake-related parameters like **--cores** and **--configfile** could be checked via 
```bash
snakemake -h
```
### Output
There are two main outputs according to the indexes used. 

**(1) target.upstream_downstream.bases.comparison.result.txt**. When the workflow used MRPP A, the main output is **target.upstream_downstream.bases.comparison.result.txt**, which contains (1) the variation locus position, (2) group1 (plus.match or minus.match, corresponding to forward-aligned reads supporting the reference allele and reverse-aligned reads supporting the reference allele), (3) group2 (plus.del or minus.del, corresponding to forward-aligned reads supporting the non-reference allele and reverse-aligned reads supporting the non-reference allele), (4) the number of reads supporting group1 (**should always be around or higher than 400 in direct sequencing**), (5) the number of reads supporting group2 (**should always be around or higher than 400 in direct sequencing**), (6) the mean current measurements of upstream and downstream config["Num"] bases centered around variation locus of group1, (7) the mean current measurements of upstream and downstream config["Num"] bases centered around variation locus of group2, (8) P values between current measurements of group1 and group2, (9) MRPP P values, (10) **MRPP A statistic, users could compare this value against the pre-set threshold (WTA/Amplicon sequencing: 0.01; direct sequencing: 0.001) in our article to decide whether the variation locus is artificial**.

**(2) fq.Qscore.info.txt**. For R9 or R10 data, when sequencing depth is low, Q score might be used to identify artificial deletions. The main output is **fq.Qscore.info.txt**, which contains (1) the variation locus position, (2) group1 (corresponding to forward-aligned reads supporting the non-reference allele and reverse-aligned reads supporting the non-reference allele), (3) group2 (plus.match or minus.match, corresponding to corresponding to forward-aligned reads supporting the reference allele and reverse-aligned reads supporting the reference allele), (4) the number of reads supporting group1 (**should always be ≥20**), (5) the number of reads supporting group2 (**should always be ≥20**), (6) the mean Q scores of upstream and downstream config["Num"] bases centered around variation locus of group1, **users could compare this value against the pre-set threshold in our article to decide whether the variation locus is artificial**, (7) the mean Q scores of upstream and downstream config["Num"] bases centered around variation locus of group2, (8) the P values between group1 and group2.

Users should note that before comparing the results to pre-set thresholds, they are strongly recommended to filter the 4th and 5th columns in the **target.upstream_downstream.bases.comparison.result.txt** or/and **fq.Qscore.info.txt** according to our article, or the identification result may be biased due to low sequencing depth.

**(3) variant.info.txt**. This file stores basic information of each variation output by VCF. The 2th-10th columns are identical to VCF. Users should note that DP4 could be lower than DP, and **choosing to use MRPP A or Q score mainly depends on DP4 field**. Our workflow use a strict criteria to extract reads with and without deletions. For example, if a deletion lacks 3 bases relative to the reference, then reads supporting the non-reference allele should only contain 3-base deletions. Therefore, some deletions may be omitted    




