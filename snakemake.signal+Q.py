configfile: "snakemake.signal+Q.config.yaml"
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


# rule bam2sam:
#     input:
#         "/public/data1/yefq/data/fast5/20220703_WGA_twist/processed/20230426_Guppy621_comparison/Sce20_guppy_sup_aligned.softclip_trimmed.endtrim10_minimap2_align.mapped.sorted.bam"
#     output:
#         "data/test/target.rev.simpified.sam.txt",
#         "data/test/target.simpified.sam.txt"
#     params:
#         opts1="--output-fmt SAM -@ 8 -h -f 16",
#         opts2="--output-fmt SAM -@ 8 -h -F 0x814",
#         region=ref_name + ":" + str(target) + "-" + str(int(target)+1)
#     shell:
#         """
#         samtools view {params.opts1} {input} {params.region} | awk -F "\\t" '{{if(NF>=10){{print $1,$2,$3,$4,$6}}}}' OFS="\\t" > {output[0]}
#         samtools view {params.opts2} {input} {params.region} | awk -F "\\t" '{{if(NF>=10){{print $1,$2,$3,$4,$6}}}}' OFS="\\t" > {output[1]}
#         """

# rule sam2readID:
#     input:
#         "data/test/target.rev.simpified.sam.txt",
#         "data/test/target.simpified.sam.txt"
#     output:
#         "data/test/target.minus.del.readID.txt",
#         "data/test/target.plus.del.readID.txt",
#         "data/test/target.minus.match.readID.txt",
#         "data/test/target.plus.match.readID.txt"
#     params:
#         opts1=str(dellen) + " " + str(target) + " D",
#         opts2=str(rglen) + " " + str(target) + " M",
#     script:
#         "scripts/Plasmid_R9-guppy_sup_Sam2ReadID.sh"

        
# rule readID2Signal:
#     input:
#         "data/test/target.minus.del.readID.txt",
#         "data/test/target.plus.del.readID.txt",
#         "data/test/target.minus.match.readID.txt",
#         "data/test/target.plus.match.readID.txt"
#     output:
#         "data/test/target.upstream_downstream.bases.singal_length.txt"
#     params:
#         opts1=str(start) + " " + str(end) + " " + str(target) + " " + str(num),
#         tombo=tombodir,
#         readIDdir="data/test"
#     script:
#         "scripts/Plasmid_R9-guppy_sup_ReadID2Signal.sh"


# rule Signal2Comparison:
#     input:
#         "data/test/target.upstream_downstream.bases.singal_length.txt"
#     output:
#         "data/test/target.upstream_downstream.bases.comparison.result.txt"
#     params:
#         subsampling=str(subsample)
#     script:
#         "scripts/Plasmid_R9-guppy_sup_Signal2Comparison.sh"






# rule Fq2info:
#     input:
#         "data/test/target.minus.del.fq",
#         "data/test/target.plus.del.fq",
#         "data/test/target.minus.match.fq",
#         "data/test/target.plus.match.fq"
#     output:
#         "data/test/target.fq.Qscore.info.txt"
#     params:
#         # Dellen=str(dellen),
#         # Target=str(target),
#         # Rglen=str(rglen)
#         opts1=str(target),
#         opts2="data/test"
#         #region="Zymo_Saccharomyces_cerevisiae_Seq5_ref:604-605"
#         #region=ref_name + ":" + str(target) + "-" + str(int(target)+1)
#     script:
#         "scripts/Plasmid_R10-guppy_sup_Fq2Summary.sh"




# rule d:
#     output:
#         done = "results/d_done",
#         md5 = "results/d_md5"
#     shell:
#         """
#         echo "done" > {output.done}

#         ssh remote -t \
#         "echo .bashrc | awk '{{cmd=\\"md5sum \\"\$0;cmd | getline md5; print md5}}' > bashrc.md5"

#         scp remote:/root/bashrc.md5 {output.md5}
#         ssh remote -t 'rm -rf bashrc.md5'
#         """









# rule all:
#     input:
#         "plots/quals.svg"


# rule bwa_map:
#     input:
#         "data/genome.fa",
#         "data/samples/{sample}.fastq"
#     output:
#         "mapped_reads/{sample}.bam"
#     shell:
#         "bwa mem {input} | samtools view -Sb - > {output}"


# rule samtools_sort:
#     input:
#         "mapped_reads/{sample}.bam"
#     output:
#         "sorted_reads/{sample}.bam"
#     shell:
#         "samtools sort -T sorted_reads/{wildcards.sample} "
#         "-O bam {input} > {output}"


# rule samtools_index:
#     input:
#         "sorted_reads/{sample}.bam"
#     output:
#         "sorted_reads/{sample}.bam.bai"
#     shell:
#         "samtools index {input}"


# rule bcftools_call:
#     input:
#         fa="data/genome.fa",
#         bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
#         bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
#     output:
#         "calls/all.vcf"
#     shell:
#         "bcftools mpileup -f {input.fa} {input.bam} | "
#         "bcftools call -mv - > {output}"


# rule plot_quals:
#     input:
#         "calls/all.vcf"
#     output:
#         "plots/quals.svg"
#     script:
#         "scripts/plot-quals.py"



# rule align:
#     input:
#         "{sample}.fq",
#         reference="ref.fa",
#     output:
#         "{sample}.sam"
#     params:
#         opts="-a -x map-ont",
#     threads: 4
#     log:
#         "align/{sample}.log"
#     conda:
#         "envs/align.yaml"  #--use-conda
#     script:
#         "scripts/align.sh"




# # #!/usr/bin/env bash

# # echo "Aligning sample ${snakemake_wildcards[sample]} with minimap2" 2> "${snakemake_log[0]}"

# # minimap2 ${snakemake_params[opts]} -t ${snakemake[threads]} "${snakemake_input[reference]}" \
# #     "${snakemake_input[0]}" > "${snakemake_output[0]}" 2>> "${snakemake_log[0]}"

