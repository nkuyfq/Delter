configfile: "R10.snakemake.config.yaml"
vcf=config["Vcf"] #240225 updated
refseq=config["Refseq"] #240225 updated
outdir=config["Outdir"] #240225 updated
BAM=config["Bam"] #240225 updated
ref_name=config["Ref"]
num=config["Num"]


print("Usage")
print("snakemake -s R10.snakemake.test.py --cores n --config Ref=refname Num=5 Vcf=path/to/VCF Refseq=path/to/refseq Outdir=path/to/outputdir Bam=path/to/sorted/bam")
print("Ref=refname".ljust(30)+"The value of #CHROM in vcf file, e.g., 'Ref=chr1'")
print("Num=5".ljust(30)+"The number of bases up- and down-stream that are centered around the variation loci, default=5")
print("Vcf=path/to/VCF".ljust(30)+"The file path to vcf file, e.g., 'Vcf=/data/res/lofreq.vcf'")
print("Refseq=path/to/refseq".ljust(30)+"The file path to reference sequence, e.g., 'Refseq=/database/COVID-19.fa'")
print("Outdir=path/to/outputdir".ljust(30)+"The file path storing the output results and intermediate files, e.g., 'Outdir=/data/res'")
print("Bam=path/to/sorted/bam".ljust(30)+"The file path to sorted bam files, e.g., 'Bam=/data/res/sorted.bam'")



rule all:
    input:
        str(outdir) + "/" + "fq.Qscore.info.txt"

# rule printhelpmessage:
#     input:
#         "data/test/helpmessage.txt"
#     output:
#         "data/test/Sce20_guppy_sup_aligned.softclip_trimmed.endtrim10_minimap2_align.mapped.sorted.bam"
#     shell:
#         """
#         cat {input}
#         """

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
        str(outdir) + "/" + "location.txt"
    script:
        "scripts/vcf2delinfo.sh"


rule delinfo2Qscoreinfo:
    input:
        str(BAM),
        str(outdir) + "/" + "targetpos.txt",
        str(outdir) + "/" + "dellen.txt",
        str(outdir) + "/" + "reglen.txt"
    output:
        str(outdir) + "/" + "fq.Qscore.info.txt"
    params:
        opts1=str(ref_name),
        opts2=str(num),
        opts3=str(outdir)
    script:
        "scripts/delinfo2Qscoreinfo.sh"


# rule bam2sam:
#     input:
#         str(BAM)
#     output:
#         "data/test/target.rev.simpified.sam.txt",
#         "data/test/target.simpified.sam.txt"
#     params:
#         opts1="--output-fmt SAM -@ 8 -h -f 16",
#         opts2="--output-fmt SAM -@ 8 -h -F 0x814",
#         #region="Zymo_Saccharomyces_cerevisiae_Seq5_ref:604-605"
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
#         "scripts/Plasmid_R10-guppy_sup_Sam2ReadID.sh"

        
# rule bam2sam_for_Fq:
#     input:
#         "data/test/Sce20_guppy_sup_aligned.softclip_trimmed.endtrim10_minimap2_align.mapped.sorted.bam"
#     output:
#         "data/test/target.rev.simpified.sam.4fq.txt",
#         "data/test/target.simpified.sam.4fq.txt"
#     params:
#         opts1="--output-fmt SAM -@ 8 -h -f 16",
#         opts2="--output-fmt SAM -@ 8 -h -F 0x814",
#         region=ref_name + ":" + str(target) + "-" + str(int(target)+1)
#     shell:
#         """
#         samtools view {params.opts1} {input} {params.region} | awk -F "\\t" '{{if(NF>=10){{print $1,$2,$4,$6,$10,$11}}}}' OFS="\\t" > {output[0]}
#         samtools view {params.opts2} {input} {params.region} | awk -F "\\t" '{{if(NF>=10){{print $1,$2,$4,$6,$10,$11}}}}' OFS="\\t" > {output[1]}
#         """

# rule readID2Fq:
#     input:
#         "data/test/target.rev.simpified.sam.4fq.txt",
#         "data/test/target.simpified.sam.4fq.txt",
#         "data/test/target.minus.del.readID.txt",
#         "data/test/target.plus.del.readID.txt",
#         "data/test/target.minus.match.readID.txt",
#         "data/test/target.plus.match.readID.txt"
#     output:
#         "data/test/target.minus.del.fq",
#         "data/test/target.plus.del.fq",
#         "data/test/target.minus.match.fq",
#         "data/test/target.plus.match.fq"
#     params:
#         opts1=str(target) + " D " + str(num),
#         opts2=str(target) + " M " + str(num),
#     script:
#         "scripts/Plasmid_R10-guppy_sup_ReadID2Fq.sh"


# rule Fq2info:
#     input:
#         "data/test/target.minus.del.fq",
#         "data/test/target.plus.del.fq",
#         "data/test/target.minus.match.fq",
#         "data/test/target.plus.match.fq"
#     output:
#         "data/test/target.fq.Qscore.info.txt"
#     params:
#         opts1=str(target),
#         opts2="data/test"
#     script:
#         "scripts/Plasmid_R10-guppy_sup_Fq2Summary.sh"



# SAMPLES = ["A", "B"]


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


