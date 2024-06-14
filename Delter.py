configfile: "Delter.config.yaml"
vcf=config["Vcf"] #240225 updated
refseq=config["Refseq"] #240225 updated
outdir=config["Outdir"] #240225 updated
BAM=config["Bam"] #240225 updated
flowcell=config["Flowcell"] #240429 updated
strategy=config["Strategy"] #240429 updated


ref_name=config["Ref"]
num=config["Num"]
tombodir=config["Tombo_dir"]
subsample=config["Subsample"]

print("Usage:")
print("snakemake -s R9.snakemake.test.py --cores 8 --config Ref=refname Num=5 Vcf=path/to/VCF Refseq=path/to/refseq Outdir=path/to/outputdir Bam=path/to/sorted/bam Tombo_dir=path/to/tombo_processed/fast5 Subsample=2000")
print("Ref=refname".ljust(40)+"The value of #CHROM in vcf file, e.g., 'Ref=chr1'")
print("Num=5".ljust(40)+"The number of bases up- and down-stream that are centered around the variation loci, default=5")
print("Vcf=path/to/VCF".ljust(40)+"The file path to vcf file, e.g., 'Vcf=/data/res/lofreq.vcf'")
print("Refseq=path/to/refseq".ljust(40)+"The file path to reference sequence, e.g., 'Refseq=/database/COVID-19.fa'")
print("Outdir=path/to/outputdir".ljust(40)+"The file path storing the output results and intermediate files, e.g., 'Outdir=/data/res'")
print("Bam=path/to/sorted/bam".ljust(40)+"The file path to sorted bam files, e.g., 'Bam=/data/res/sorted.bam'")
print("Tombo_dir=path/to/tombo_processed/fast5".ljust(40)+"The file path to tombo-resquiggled single fats5 files, e.g., 'Tombo_dir=/data/fast5'")
print("Subsample=2000".ljust(40)+"The number to subsample from reads covering variation loci, should be larger than 200, default=2000")
print("Flowcell=R9".ljust(40)+"The version of flow cell, should be R9 or R10, default=R9")
print("Strategy=Amplicon".ljust(40)+"The sequencing strategy, should be Amplicon or Direct, default=Amplicon")



rule all:
    input:
        str(outdir) + "/" + "run.log"

rule vcf2delinfo:
    input:
        str(refseq),
        str(vcf)
    output:
        str(outdir) + "/" + "variant.info.txt",
        str(outdir) + "/" + "targetpos.txt",
        str(outdir) + "/" + "startpos.txt",
        str(outdir) + "/" + "endpos.txt",
        str(outdir) + "/" + "dellen.txt",
        str(outdir) + "/" + "reglen.txt",
        str(outdir) + "/" + "location.txt",
        str(outdir) + "/" + "metric.txt"
    params:
        opts1=str(flowcell),
        opts2=str(strategy)
    script:
        "scripts/vcf2delinfo-v2.sh"

rule delinfo2Signal_Qinfo:
    input:
        str(BAM),
        str(outdir) + "/" + "targetpos.txt",
        str(outdir) + "/" + "startpos.txt",
        str(outdir) + "/" + "endpos.txt",
        str(outdir) + "/" + "dellen.txt",
        str(outdir) + "/" + "reglen.txt",
        str(outdir) + "/" + "metric.txt"
    output:
        str(outdir) + "/" + "target.upstream_downstream.bases.comparison.result.txt",
        str(outdir) + "/" + "fq.Qscore.info.txt",
        str(outdir) + "/" + "run.log"
    params:
        opts1=str(ref_name),
        opts2=str(num),
        opts3=str(outdir),
        opts4=str(tombodir),
        opts5=str(subsample)
    script:
        "scripts/delinfo2Signal+Qinfo.sh"

