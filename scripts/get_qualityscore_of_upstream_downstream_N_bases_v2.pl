#!/usr/bin/perl -w
use warnings;
use strict;


#updated on 20230715 
#以ReadID.txt为参考，在ReadID.txt的reads导出相应fastq，不在txt里的，不导出。


#sam文件中的$1,$2,$4,$6,$10,$11:queryname,flag,position,CIGAR,queryseq,qualityscore（如果比对到负链，queryseq为原始read的反向互补序列，qualityscore为反向输出）


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


my $targetpos=$ARGV[0];
my $cigarOP=$ARGV[1];
my $substrlen=$ARGV[2];
my $samfile=$ARGV[3];
my $output=$ARGV[4];
my $listfile=$ARGV[5];

#sam文件中的$1,$2,$4,$6,$10,$11:queryname,flag,position,CIGAR,queryseq,qualityscore（如果比对到负链，queryseq为原始read的反向互补序列，qualityscore为反向输出）

my @readids;
open FP,"$listfile";
while(<FP>){
	chomp;
	push @readids,$_;
}
close FP;

open FP,"< $samfile";
open FP1, ">$output";
while(<FP>){
	chomp;
	my @tmps=split/\t/;
	if(is_ele($tmps[0],\@readids)){
		my $seqname=$tmps[0];
		my $flag=$tmps[1];
		my $startpos=$tmps[2];
		my $cigar=$tmps[3];
		my $query=$tmps[4];
		my $quality=$tmps[5];
		my $curcoord=$startpos-1; #coord in reference, 1-based
		my $querycoord=0;	  #coord in query sequence, 1-based
		my @matchs = ($cigar=~ /(\d+)([SMDI]+)/g);
		my @nums=($cigar=~ /(\d+)/g);
		my @ops=($cigar=~ /([SMDI]+)/g);
		my $substrstartpos=0;	  #substring start position, 1-based
		for(my $i=0;$i<@nums;$i++){
			if($ops[$i] eq "S"){
				$querycoord=$querycoord+$nums[$i];
			}
			if($ops[$i] eq "I"){
				$querycoord=$querycoord+$nums[$i];
			}
			if($ops[$i] eq "D"){
				$curcoord=$curcoord+$nums[$i];
			}
			if($ops[$i] eq "M"){
				$curcoord=$curcoord+$nums[$i];
				$querycoord=$querycoord+$nums[$i];
			}
			if($curcoord < $targetpos){
			}
			if($curcoord == $targetpos){
				if($ops[$i] eq $cigarOP){
					$substrstartpos=$querycoord;
					last;
				}
				else {
					last;
				}
			}
			if($curcoord > $targetpos){
				if($ops[$i] eq $cigarOP){
					if($ops[$i] eq "D"){
						$substrstartpos=$querycoord;
						last;
					}
					if($ops[$i] eq "M"){
						$substrstartpos=$querycoord-($curcoord-$targetpos);
						last;
					}
				}
				else {
					last;
				}
			}

		}
		if(($substrstartpos-$substrlen >= 1) && ($substrstartpos+$substrlen <= length($query))){
			my $idx=$substrstartpos-$substrlen-1;
			my $len=2*$substrlen+1;
			if($flag == 16){
				my $subseq=substr($query,$idx,$len);
				my $subquality=substr($quality,$idx,$len);
				$subquality=reverse($subquality);
				$subseq=reverse_complement($subseq);
				print FP1 "@".$seqname."\n$subseq\n+\n$subquality\n";
			}
			if($flag == 0){
				my $subseq=substr($query,$idx,$len);
				my $subquality=substr($quality,$idx,$len);
				print FP1 "@".$seqname."\n$subseq\n+\n$subquality\n";
			}
		}
	}
	
}
close FP;
close FP1;


print "Job is done\n";
