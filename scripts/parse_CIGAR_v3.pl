#!/usr/bin/perl -w
use warnings;
use strict;

#updated on 20230713,针对CIGAR为“M”时的情况设计，
#将Rglen纳入考虑范围，比如pos123处的Rglen为3，那么只输出在pos123-125处都为"M"的ReadID
#perl /data1/yefq/scripts/yfq_remote/parse_CIGAR.pl 9812 D in.sam readID.txt


sub is_ele{
	my $factor=$_[0];
	my $arr=$_[1];
	my $flag=0;
	foreach my $arr_mem (@$arr){
		if ($factor eq $arr_mem){
			$flag=1;
			last;
		}
	}
	$flag;
}


sub reverse_complement {
        my $dna = $_[0];
        # reverse the DNA sequence
        my $revcomp = reverse($dna);
        # complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}


my $Rglen=$ARGV[0];
my $refpos=$ARGV[1];
my $cigarOP=$ARGV[2];
my $samfile=$ARGV[3];
my $output=$ARGV[4];


my @readids;
open FP,"< $samfile";
while(<FP>){
	chomp;
	my @tmps=split/\t/;
	my $seqname=$tmps[0];
	my $flag=$tmps[1];
	my $startpos=$tmps[3];
	my $cigar=$tmps[4];
	my $curcoord=$startpos-1;
	my @matchs = ($cigar=~ /(\d+)([SMDI]+)/g);
	my @nums=($cigar=~ /(\d+)/g);
	my @ops=($cigar=~ /([SMDI]+)/g);
	for(my $i=0;$i<@nums;$i++){
		if(($ops[$i] ne "S") && ($ops[$i] ne "I")){
			$curcoord=$curcoord+$nums[$i];
			if($curcoord < $refpos){
			}
			if($curcoord == $refpos){
				last
				# if($ops[$i] eq $cigarOP){   #此时，下一个CIGAR可能是D或者I，如果下一个是D，这不考虑此ReadID；如果下一个是I，下下个是M，且$nums[$i]>=$Rglen-1，则考虑此ReadID
				# 	if($i+2 < @nums && $ops[$i+1] == "I" && $ops[$i+2] == "M" && $nums[$i+2] >= ($Rglen-1)){
				# 		print "Catched\n";
				# 		push @readids,$seqname;
				# 		last;
				# 	}	
				# }
				# else {
				# 	last;
				# }
			}
			if($curcoord > $refpos){
				if($ops[$i] eq $cigarOP && $curcoord-$refpos >= ($Rglen-1)){
					#print "Hey, Catched\n";
					push @readids,$seqname;
					last;
				}
				else {
					last;
				}
			}
		}
	}
}
close FP;

open FP,"> $output";
foreach my $id (@readids){
	print FP "$id\n";
}
close FP;



print "Job is done\n";
