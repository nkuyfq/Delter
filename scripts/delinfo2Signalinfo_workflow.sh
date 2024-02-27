#!/usr/bin/env bash
export HDF5_PLUGIN_PATH=lib/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin

bam=$1
targetposlist=$2
startposlist=$3
endposlist=$4
dellenlist=$5
reglenlist=$6
refname=$7
num=$8
outdir=$9
tombodir=${10}
subsample=${11}
output=${12}

targetposarray1=()
for line in $(cat $targetposlist)
do
	targetposarray1+=($line)
	#echo $line
done

startposarray1=()
for line in $(cat $startposlist)
do
	startposarray1+=($line)
	#echo $line
done

endposarray1=()
for line in $(cat $endposlist)
do
	endposarray1+=($line)
	#echo $line
done

Dellens=()
for line in $(cat $dellenlist)
do
	Dellens+=($line)
	#echo $line
done

Rglens=()
for line in $(cat $reglenlist)
do
	Rglens+=($line)
	#echo $line
done

arraylen=${#targetposarray1[*]}

#====output readIDs====
echo "Output ReadIDs"
for ((i=0; i<$arraylen; i ++))
do
	targetpos=${targetposarray1[$i]}
	Dellen=${Dellens[$i]}
	Rglen=${Rglens[$i]}
	echo ${targetpos}
	nextpos=$(($targetpos+1))
	samtools view --output-fmt SAM -@ 8 -h -f 16 ${bam} $refname:${targetpos}-${nextpos} | awk -F "\t" '{if(NF>=10){print $1,$2,$3,$4,$6}}' OFS="\t" > ${outdir}/pos${targetpos}.rev.simpified.sam.txt
	samtools view --output-fmt SAM -@ 8 -h -F 0x814 ${bam} $refname:${targetpos}-${nextpos} | awk -F "\t" '{if(NF>=10){print $1,$2,$3,$4,$6}}' OFS="\t" > ${outdir}/pos${targetpos}.simpified.sam.txt
	perl scripts/parse_CIGAR_v2.pl ${Dellen} ${targetpos} D ${outdir}/pos${targetpos}.rev.simpified.sam.txt ${outdir}/target.pos${targetpos}.minus.del.readID.txt
	perl scripts/parse_CIGAR_v2.pl ${Dellen} ${targetpos} D ${outdir}/pos${targetpos}.simpified.sam.txt ${outdir}/target.pos${targetpos}.plus.del.readID.txt
	perl scripts/parse_CIGAR_v3.pl ${Rglen} ${targetpos} M ${outdir}/pos${targetpos}.rev.simpified.sam.txt ${outdir}/target.pos${targetpos}.minus.match.readID.txt
	perl scripts/parse_CIGAR_v3.pl ${Rglen} ${targetpos} M ${outdir}/pos${targetpos}.simpified.sam.txt ${outdir}/target.pos${targetpos}.plus.match.readID.txt
	rm ${outdir}/*.simpified.sam.txt
done

#====output signal====
echo "Output signal"
for ((i=0; i<$arraylen; i ++))
do
	targetpos=${targetposarray1[$i]}
	startpos=${startposarray1[$i]}
	endpos=${endposarray1[$i]}
	echo ${targetpos}
	#using a new version of python script
	python scripts/get_signals_of_unpstream_downstream_N_bases_v3-240226.py ${startpos} ${endpos} ${targetpos} ${num} ${tombodir} ${outdir} ${outdir} 
done

#====compare signal====
echo "Comparing signal"

Rscript scripts/comparison_Del_upstream_downstream_N_bases_subsample_v4-240226.R ${subsample} ${outdir} ${output}



# 	perl $zilinplace/scripts/yfq_remote/get_qualityscore_of_upstream_downstream_N_bases_v2.pl ${targetpos} D $num ${name}_guppy_sup.pos${targetpos}.rev.simpified.sam.txt ./Fq/TP/${name}_guppy_sup.pos${targetpos}.minus.del.fq ${name}_guppy_sup.pos${targetpos}.minus.del.readID.txt
# 	perl $zilinplace/scripts/yfq_remote/get_qualityscore_of_upstream_downstream_N_bases_v2.pl ${targetpos} D $num ${name}_guppy_sup.pos${targetpos}.simpified.sam.txt ./Fq/TP/${name}_guppy_sup.pos${targetpos}.plus.del.fq ${name}_guppy_sup.pos${targetpos}.plus.del.readID.txt
# 	perl $zilinplace/scripts/yfq_remote/get_qualityscore_of_upstream_downstream_N_bases_v2.pl ${targetpos} M $num ${name}_guppy_sup.pos${targetpos}.rev.simpified.sam.txt ./Fq/TP/${name}_guppy_sup.pos${targetpos}.minus.match.fq ${name}_guppy_sup.pos${targetpos}.minus.match.readID.txt
# 	perl $zilinplace/scripts/yfq_remote/get_qualityscore_of_upstream_downstream_N_bases_v2.pl ${targetpos} M $num ${name}_guppy_sup.pos${targetpos}.simpified.sam.txt ./Fq/TP/${name}_guppy_sup.pos${targetpos}.plus.match.fq ${name}_guppy_sup.pos${targetpos}.plus.match.readID.txt


# # minimap2 ${snakemake_params[opts]} -t ${snakemake[threads]} "${snakemake_input[reference]}" \
# #     "${snakemake_input[0]}" > "${snakemake_output[0]}" 2>> "${snakemake_log[0]}"


# #==============Step 5.11 Sce20: parse bam files to get readIDs and output current signals==================
# name="Sce20"
# #first parse bam files to get readIDs
# targetposarray1=("604" "661" "685" "85" "250" "354" "460" "286" "327" "816" "625" "1523" "981" "1366")
# Dellens=("1" "1" "1" "1" "1" "2" "2" "1" "1" "2" "2" "3" "1" "1")
# Rglens=("4" "5" "4" "1" "3" "2" "2" "1" "1" "2" "7" "12" "5" "4")

# arraylen=${#targetposarray1[*]}

# conda activate evaluation
# for ((i=0; i<$arraylen; i ++))
# do
# 	targetpos=${targetposarray1[$i]}
# 	Dellen=${Dellens[$i]}
# 	Rglen=${Rglens[$i]}
# 	echo ${targetpos}
# 	nextpos=$(($targetpos+1))
# 	cd $rawdata/data_process/guppy_sup
# 	samtools view --output-fmt SAM -@ 8 -h -f 16 ${name}_guppy_sup_aligned.softclip_trimmed.endtrim10_minimap2_align.mapped.sorted.bam Zymo_Saccharomyces_cerevisiae_Seq5_ref:${targetpos}-${nextpos} | awk -F "\t" '{if(NF>=10){print $1,$2,$3,$4,$6}}' OFS="\t" > ${name}_guppy_sup.pos${targetpos}.rev.simpified.sam.txt
# 	samtools view --output-fmt SAM -@ 8 -h -F 0x814 ${name}_guppy_sup_aligned.softclip_trimmed.endtrim10_minimap2_align.mapped.sorted.bam Zymo_Saccharomyces_cerevisiae_Seq5_ref:${targetpos}-${nextpos} | awk -F "\t" '{if(NF>=10){print $1,$2,$3,$4,$6}}' OFS="\t" > ${name}_guppy_sup.pos${targetpos}.simpified.sam.txt
# 	perl $zilinplace/scripts/yfq_remote/parse_CIGAR_v2.pl ${Dellen} ${targetpos} D ${name}_guppy_sup.pos${targetpos}.rev.simpified.sam.txt ${name}_guppy_sup.pos${targetpos}.minus.del.readID.txt
# 	perl $zilinplace/scripts/yfq_remote/parse_CIGAR_v2.pl ${Dellen} ${targetpos} D ${name}_guppy_sup.pos${targetpos}.simpified.sam.txt ${name}_guppy_sup.pos${targetpos}.plus.del.readID.txt
# #	perl $zilinplace/scripts/yfq_remote/parse_CIGAR.pl ${targetpos} D ${name}_guppy_sup.pos${targetpos}.rev.simpified.sam.txt ${name}_guppy_sup.pos${targetpos}.minus.del.readID.txt
# #	perl $zilinplace/scripts/yfq_remote/parse_CIGAR.pl ${targetpos} D ${name}_guppy_sup.pos${targetpos}.simpified.sam.txt ${name}_guppy_sup.pos${targetpos}.plus.del.readID.txt
# 	perl $zilinplace/scripts/yfq_remote/parse_CIGAR_v3.pl ${Rglen} ${targetpos} M ${name}_guppy_sup.pos${targetpos}.rev.simpified.sam.txt ${name}_guppy_sup.pos${targetpos}.minus.match.readID.txt
# 	perl $zilinplace/scripts/yfq_remote/parse_CIGAR_v3.pl ${Rglen} ${targetpos} M ${name}_guppy_sup.pos${targetpos}.simpified.sam.txt ${name}_guppy_sup.pos${targetpos}.plus.match.readID.txt
# done


# #===============================Step 6 output fastq of FP/TP positions=================================

# #==============Step 6.11 Sce20: parse bam files to get readIDs and output current signals==================
# name="Sce20"
# #first parse bam files to get readIDs
# conda activate evaluation
# for targetpos in {"604","661","685","85","250","354","460","286","327","816","625","1523"}
# do
# 	echo ${targetpos}
# 	nextpos=$(($targetpos+1))
# 	cd $rawdata/data_process/guppy_sup
# 	samtools view --output-fmt SAM -@ 8 -h -f 16 ${name}_guppy_sup_aligned.softclip_trimmed.endtrim10_minimap2_align.mapped.sorted.bam Zymo_Saccharomyces_cerevisiae_Seq5_ref:${targetpos}-${nextpos} | awk -F "\t" '{if(NF>=10){print $1,$2,$4,$6,$10,$11}}' OFS="\t" > ${name}_guppy_sup.pos${targetpos}.rev.simpified.sam.txt
# 	samtools view --output-fmt SAM -@ 8 -h -F 0x814 ${name}_guppy_sup_aligned.softclip_trimmed.endtrim10_minimap2_align.mapped.sorted.bam Zymo_Saccharomyces_cerevisiae_Seq5_ref:${targetpos}-${nextpos} | awk -F "\t" '{if(NF>=10){print $1,$2,$4,$6,$10,$11}}' OFS="\t" > ${name}_guppy_sup.pos${targetpos}.simpified.sam.txt
# 	perl $zilinplace/scripts/yfq_remote/get_qualityscore_of_upstream_downstream_N_bases_v2.pl ${targetpos} D $num ${name}_guppy_sup.pos${targetpos}.rev.simpified.sam.txt ./Fq/TP/${name}_guppy_sup.pos${targetpos}.minus.del.fq ${name}_guppy_sup.pos${targetpos}.minus.del.readID.txt
# 	perl $zilinplace/scripts/yfq_remote/get_qualityscore_of_upstream_downstream_N_bases_v2.pl ${targetpos} D $num ${name}_guppy_sup.pos${targetpos}.simpified.sam.txt ./Fq/TP/${name}_guppy_sup.pos${targetpos}.plus.del.fq ${name}_guppy_sup.pos${targetpos}.plus.del.readID.txt
# 	perl $zilinplace/scripts/yfq_remote/get_qualityscore_of_upstream_downstream_N_bases_v2.pl ${targetpos} M $num ${name}_guppy_sup.pos${targetpos}.rev.simpified.sam.txt ./Fq/TP/${name}_guppy_sup.pos${targetpos}.minus.match.fq ${name}_guppy_sup.pos${targetpos}.minus.match.readID.txt
# 	perl $zilinplace/scripts/yfq_remote/get_qualityscore_of_upstream_downstream_N_bases_v2.pl ${targetpos} M $num ${name}_guppy_sup.pos${targetpos}.simpified.sam.txt ./Fq/TP/${name}_guppy_sup.pos${targetpos}.plus.match.fq ${name}_guppy_sup.pos${targetpos}.plus.match.readID.txt
# 	#cat ${name}_guppy_sup.pos${targetpos}.*.fq > ${name}_guppy_sup.pos${targetpos}.combined.fq
# done
# conda deactivate
# #cat ${name}_guppy_sup.pos*.combined.fq > ${name}_guppy_sup.allpos.combined.fq
