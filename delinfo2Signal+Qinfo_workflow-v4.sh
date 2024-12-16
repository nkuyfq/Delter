#!/usr/bin/env bash
export HDF5_PLUGIN_PATH=lib/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin

bam=$1
targetposlist=$2
startposlist=$3
endposlist=$4
dellenlist=$5
reglenlist=$6
metriclist=$7
refname=$8
num=$9
outdir=${10}
tombodir=${11}
subsample=${12}
output1=${13}
output2=${14}
output3=${15}
usermetric=${16}


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

Metrics=()
for line in $(cat $metriclist)
do
	Metrics+=($line)
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

#====output signal or fastqs====

declare -i signalcount=0
declare -i qscorecount=0

for ((i=0; i<$arraylen; i ++))
do
	targetpos=${targetposarray1[$i]}
	startpos=${startposarray1[$i]}
	endpos=${endposarray1[$i]}
	metric=${Metrics[$i]}
	echo ${targetpos}
	nextpos=$(($targetpos+1))
	if [ $usermetric = "NA" ]; then
		if [ $metric = "Signal" ]; then
			#using a new version of python script
			echo "Extracting sequencing current signals"
			signalcount=signalcount+1
			python scripts/get_signals_of_unpstream_downstream_N_bases_v3-240226.py ${startpos} ${endpos} ${targetpos} ${num} ${tombodir} ${outdir} ${outdir}
		fi
		if [ $metric = "Qscore" ]; then
			echo "Extracting base qualities"
			qscorecount=qscorecount+1
			samtools view --output-fmt SAM -@ 8 -h -f 16 ${bam} $refname:${targetpos}-${nextpos} | awk -F "\t" '{if(NF>=10){print $1,$2,$4,$6,$10,$11}}' OFS="\t" > ${outdir}/pos${targetpos}.rev.simpified.sam.txt
			samtools view --output-fmt SAM -@ 8 -h -F 0x814 ${bam} $refname:${targetpos}-${nextpos} | awk -F "\t" '{if(NF>=10){print $1,$2,$4,$6,$10,$11}}' OFS="\t" > ${outdir}/pos${targetpos}.simpified.sam.txt
			perl scripts/get_qualityscore_of_upstream_downstream_N_bases_v2.pl ${targetpos} D $num ${outdir}/pos${targetpos}.rev.simpified.sam.txt ${outdir}/target.pos${targetpos}.minus.del.fq ${outdir}/target.pos${targetpos}.minus.del.readID.txt
			perl scripts/get_qualityscore_of_upstream_downstream_N_bases_v2.pl ${targetpos} D $num ${outdir}/pos${targetpos}.simpified.sam.txt ${outdir}/target.pos${targetpos}.plus.del.fq ${outdir}/target.pos${targetpos}.plus.del.readID.txt
			perl scripts/get_qualityscore_of_upstream_downstream_N_bases_v2.pl ${targetpos} M $num ${outdir}/pos${targetpos}.rev.simpified.sam.txt ${outdir}/target.pos${targetpos}.minus.match.fq ${outdir}/target.pos${targetpos}.minus.match.readID.txt
			perl scripts/get_qualityscore_of_upstream_downstream_N_bases_v2.pl ${targetpos} M $num ${outdir}/pos${targetpos}.simpified.sam.txt ${outdir}/target.pos${targetpos}.plus.match.fq ${outdir}/target.pos${targetpos}.plus.match.readID.txt
			if [ -s ${outdir}/target.pos${targetpos}.minus.del.fq -a -s ${outdir}/target.pos${targetpos}.minus.match.fq -a -s ${outdir}/target.pos${targetpos}.plus.del.fq -a -s ${outdir}/target.pos${targetpos}.plus.match.fq ]; then
				echo "Files OK"
			else
				rm ${outdir}/target.pos${targetpos}.minus.del.fq
				rm ${outdir}/target.pos${targetpos}.minus.match.fq
				rm ${outdir}/target.pos${targetpos}.plus.del.fq
				rm ${outdir}/target.pos${targetpos}.plus.match.fq
			fi
			rm ${outdir}/*.simpified.sam.txt
		fi
	elif [ $usermetric = "Signal" ]; then
		echo "Extracting sequencing current signals"
		signalcount=signalcount+1
		python scripts/get_signals_of_unpstream_downstream_N_bases_v3-240226.py ${startpos} ${endpos} ${targetpos} ${num} ${tombodir} ${outdir} ${outdir}
	elif [ $usermetric = "Qscore" ]; then
		echo "Extracting base qualities"
		qscorecount=qscorecount+1
		samtools view --output-fmt SAM -@ 8 -h -f 16 ${bam} $refname:${targetpos}-${nextpos} | awk -F "\t" '{if(NF>=10){print $1,$2,$4,$6,$10,$11}}' OFS="\t" > ${outdir}/pos${targetpos}.rev.simpified.sam.txt
		samtools view --output-fmt SAM -@ 8 -h -F 0x814 ${bam} $refname:${targetpos}-${nextpos} | awk -F "\t" '{if(NF>=10){print $1,$2,$4,$6,$10,$11}}' OFS="\t" > ${outdir}/pos${targetpos}.simpified.sam.txt
		perl scripts/get_qualityscore_of_upstream_downstream_N_bases_v2.pl ${targetpos} D $num ${outdir}/pos${targetpos}.rev.simpified.sam.txt ${outdir}/target.pos${targetpos}.minus.del.fq ${outdir}/target.pos${targetpos}.minus.del.readID.txt
		perl scripts/get_qualityscore_of_upstream_downstream_N_bases_v2.pl ${targetpos} D $num ${outdir}/pos${targetpos}.simpified.sam.txt ${outdir}/target.pos${targetpos}.plus.del.fq ${outdir}/target.pos${targetpos}.plus.del.readID.txt
		perl scripts/get_qualityscore_of_upstream_downstream_N_bases_v2.pl ${targetpos} M $num ${outdir}/pos${targetpos}.rev.simpified.sam.txt ${outdir}/target.pos${targetpos}.minus.match.fq ${outdir}/target.pos${targetpos}.minus.match.readID.txt
		perl scripts/get_qualityscore_of_upstream_downstream_N_bases_v2.pl ${targetpos} M $num ${outdir}/pos${targetpos}.simpified.sam.txt ${outdir}/target.pos${targetpos}.plus.match.fq ${outdir}/target.pos${targetpos}.plus.match.readID.txt
		if [ -s ${outdir}/target.pos${targetpos}.minus.del.fq -a -s ${outdir}/target.pos${targetpos}.minus.match.fq -a -s ${outdir}/target.pos${targetpos}.plus.del.fq -a -s ${outdir}/target.pos${targetpos}.plus.match.fq ]; then
			echo "Files OK"
		else
			rm ${outdir}/target.pos${targetpos}.minus.del.fq
			rm ${outdir}/target.pos${targetpos}.minus.match.fq
			rm ${outdir}/target.pos${targetpos}.plus.del.fq
			rm ${outdir}/target.pos${targetpos}.plus.match.fq
		fi
		rm ${outdir}/*.simpified.sam.txt
	else
		echo "Error: User provided metric is wrong, please check!"
	fi

done

#====compare signal or fastq====
if [ $signalcount -ge 1 ]; then
	echo "Comparing signal"
	Rscript scripts/comparison_Del_upstream_downstream_N_bases_subsample_v4-240226.R ${subsample} ${outdir} ${output1}
	ls ${outdir} | grep "singal_length.txt" | xargs rm 
else
	echo "No need to compare current signals" > ${output1}
fi

if [ $qscorecount -ge 1 ]; then
	echo "Comparing fastq"
	Rscript scripts/R10_TP_FP_Qscore.R ${outdir} ${output2}
	ls ${outdir} | grep "minus.del.fq" | xargs rm
	ls ${outdir} | grep "minus.match.fq" | xargs rm
	ls ${outdir} | grep "plus.del.fq" | xargs rm
	ls ${outdir} | grep "plus.match.fq" | xargs rm
else
	echo "No need to compare fastqs" > ${output2}
fi

echo "Run completed" > ${output3}


ls ${outdir} | grep "readID.txt" | xargs rm 
rm $targetposlist
rm $startposlist
rm $endposlist
rm $dellenlist
rm $reglenlist
#rm $metriclist

