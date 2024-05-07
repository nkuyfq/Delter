#!/usr/bin/perl -w
use warnings;
use strict;

#perl parse_variant_type_via_lofreq_output.pl Zymo_Saccharomyces_cerevisiae_Seq5_ref.fa Sce20_guppy_sup_aligned.softclip_trimmed.endtrim10_minimap2_align.lofreq.filtered.processed.vcf variant.info.txt targetpos.txt startpos.txt endpos.txt dellen.txt reglen.txt location.txt


my $refseq=$ARGV[0];
my $input=$ARGV[1]; #vcf file
my $flowcell=$ARGV[2];
my $strategy=$ARGV[3];
my $output=$ARGV[4];
my $targetposfile=$ARGV[5];
my $startposfile=$ARGV[6];
my $endposfile=$ARGV[7];
my $dellenfile=$ARGV[8];
my $reglenfile=$ARGV[9];
my $locationfile=$ARGV[10];
my $metricfile=$ARGV[11];


#define coverage threshold
my %threshold;
$threshold{"R9_Amplicon"}=20;
$threshold{"R9_Direct"}=400;
$threshold{"R10_Direct"}=20;
$threshold{"R10_Amplicon"}=20;


open FP, "< $refseq";
my @ref;
while(<FP>){
	chomp;
	if(/^>(.*)/){

	}
	else {
		my @tmp=split//,$_;
		push @ref,@tmp;
	}
}
close FP;
print "Reference is loaded.\n";
#print "$ref[0]\n";

open FP2, ">$output";
print FP2 "Variant\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tDP\tAF\tDP4\tType\tLocation\tTargetPos\tSTARTPOS\tENDPOS\tDellen\tRegionlen\n";

open FP3, ">$targetposfile";
open FP4, ">$startposfile";
open FP5, ">$endposfile";
open FP6, ">$dellenfile";
open FP7, ">$reglenfile";
open FP8, ">$locationfile";
open FP9, ">$metricfile";


#print "$curfile\n";
# my @tmps=split/\//,$curfile;
# my $filename=$tmps[@tmps-1];
# my $variant=$kword;
my $c=0;
open FP1, "<$input";
while(<FP1>){
	chomp;
	$c++;
	if($c % 10 eq 0){
		print "$c\n";
	}
	my $POS;
	my $DP;
	my $AF;
	my $DP4;
	my $Type;
	my $Location;
	my $variant;
	my $TargetPos;
	my $STARTPOS;
	my $ENDPOS;
	my $Regionlen;
	my $Dellen;
	my $Metric;
	if(!/^#/){
		my @infos=split/\t/,$_;
		$variant=$infos[0];
		#my @ref=@{$refs{$variant}};
		$POS=$infos[1];  #0-based
		my $refbase=$infos[3];
		my $altbase=$infos[4];
		if(length($refbase) < length($altbase)){
			$Type="INDEL";
			$Location="Ins";
			if($infos[7] =~ /DP=(.*);AF=(.*);SB=(.*);DP4=(.*);INDEL/){
				$DP=$1;
				$AF=$2;
				$DP4=$4;
				$TargetPos=0;
				$STARTPOS=0;
				$ENDPOS=0;
				$Regionlen=0;
				$Dellen=0;
				$Metric=0;
			}
		}
		if(length($refbase) == length($altbase)){
			$Type="SNP";
			$Location="SNP";
			if($infos[7] =~ /DP=(.*);AF=(.*);SB=(.*);DP4=(.*)/){
				$DP=$1;
				$AF=$2;
				$DP4=$4;
				$TargetPos=0;
				$STARTPOS=0;
				$ENDPOS=0;
				$Regionlen=0;
				$Dellen=0;
				$Metric=0;
			}
		}
		if(length($refbase) > length($altbase)){
			$Type="INDEL";
			my $i=$POS;
			while($ref[$POS] eq $ref[$i]){
				$i++;
			}
			if($i-$POS ge 3){
				$Location="Homo";
				$TargetPos=$POS+1;
				$STARTPOS=$POS;
				$ENDPOS=$i-1;
				$Regionlen=$i-$POS;
				$Dellen=length($refbase)-length($altbase);

			}
			else {
				$Location="Non-Homo";
				$TargetPos=$POS+1;
				$STARTPOS=$POS;
				$ENDPOS=$POS+length($refbase)-length($altbase)-1;
				$Regionlen=length($refbase)-length($altbase);
				$Dellen=length($refbase)-length($altbase);
			}
			if($infos[7] =~ /DP=(.*);AF=(.*);SB=(.*);DP4=(.*);INDEL/){
				$DP=$1;
				$AF=$2;
				$DP4=$4;
				if(($flowcell eq "R9") && ($strategy eq "Amplicon")){
					$Metric="Signal";
				}
				if($flowcell eq "R10"){
					$Metric="Qscore";
				}
				if(($flowcell eq "R9") && ($strategy eq "Direct")){
					my $str=$flowcell."_".$strategy;
					my $curthres=$threshold{$str};
					my @tmps=split/,/,$DP4;
					if($tmps[2] < $curthres && $tmps[3] < $curthres){
						$Metric="Qscore";
					}
					else {
						$Metric="Signal";
					}
				}
			}
			print FP2 "$variant\t$infos[1]\t$infos[2]\t$infos[3]\t$infos[4]\t$infos[5]\t$infos[6]\t$DP\t$AF\t$DP4\t$Type\t$Location\t$TargetPos\t$STARTPOS\t$ENDPOS\t$Dellen\t$Regionlen\n";
			print FP3 "$TargetPos\n";
			print FP4 "$STARTPOS\n";
			print FP5 "$ENDPOS\n";
			print FP6 "$Dellen\n";
			print FP7 "$Regionlen\n";
			print FP8 "$Location\n";
			print FP9 "$Metric\n";
# open FP3, ">$targetposfile";
# open FP4, ">$startposfile";
# open FP5, ">$endposfile";
# open FP6, ">$dellenfile";
# open FP7, ">$reglenfile";

		}
	}
}
close FP1;
close FP2;
close FP3;
close FP4;
close FP5;
close FP6;
close FP7;
close FP8;
close FP9;

print "Job is done.\n";




# open FP, "< $file" or die "$!\n";
# open FP1, "> $output" or die "$!\n";
# my $seq_num=0;
# my $header_flag=0;
# my $cur_length=0;
# my @headers;
# #my @taxnames;
# my @lengths;

# while(<FP>){
# 	chomp;	
# 	if(/^>(.*)/){
# 		$seq_num++;
# 		$header_flag=1;
# 		#print "seq_num is $seq_num\n";
# 		push @headers,$1;
# 		push @lengths,$cur_length;
# 		$cur_length=0;
# 	}
# 	else {
# 		$header_flag=0;
# 	}
	
# 	if($header_flag == 0){
# 		$cur_length=$cur_length+length($_);
# 	}
# }
# push @lengths,$cur_length;

# close FP;

# my @real_len=@lengths[1..$#lengths];

# my $num=@headers;
# #my $num1=@taxnames;
# my $num2=@lengths;
# my $num3=@real_len;

# #print "seq number  is $num\t $num1\t $num2\t$num3\n";

# for(my $i=0;$i<=$num-1;$i++){
# 	print FP1 "$headers[$i]\t$real_len[$i]\n";
# }
# close FP1;

# print "Job is done";


