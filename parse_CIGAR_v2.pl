#!/usr/bin/perl -w
use warnings;
use strict;

#updated on 20230712
#将Dellen纳入考虑范围，比如pos123处的Dellen为3，那么只输出在pos123-125处都有Del，且在pos126处无Del的ReadID
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


my $Dellen=$ARGV[0];
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
				if($ops[$i] eq $cigarOP && $Dellen == 1 && $nums[$i] == 1){
					#print "Catched\n";
					push @readids,$seqname;
					last;
				}
				else {
					last;
				}
			}
			if($curcoord > $refpos){
				if($ops[$i] eq $cigarOP && $nums[$i] == $Dellen && $curcoord-$refpos == ($Dellen-1)){
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
