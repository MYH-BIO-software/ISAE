#!/usr/bin/perl -w
my $begin_time=time;
use strict;
use Getopt::Long;
my ($output,$fr,$rev);
GetOptions(
			"help|?" =>\&USAGE, 
			"out:s"=>\$output,
			"fr:s"=>\$fr,
			"rev:s"=>\$rev
				) or &USAGE;
&USAGE unless ($output && $fr && $rev) ;
sub USAGE{
	my $usage=<<"USAGE";
	Description:
	This program can merge wellington footprint analysis wig files.
	Usage:
	-fr				input file				Must be given, the fr wig file.
	-rev				input file				Must be given, the rev wig file.
	-out				output file path				output file path must be given. 
USAGE
	print $usage;
	exit;
}
#mkdir $output unless (-d $output); 
my @perl=split /\//,$0; my $sh_output=$output.".sh";
open SH, ">".$sh_output or die $sh_output,"\t$!";
print SH join(" ","perl",$0,"-fr",$fr,"-rev",$rev,"-out",$output)," ";
close SH;
####
my $tmp="$output\.tmp";
`paste $fr $rev > $tmp`;
open IN, $tmp or die $tmp," $!\n";
open OUT, ">".$output or die $output,"\t$!\n";
open OUT2, ">".$rev.".plus.wig" or die $output,"\t$!\n";
while(<IN>){
	chomp;
	$_=~s/start\=0/start\=1/ig if($_=~/[A-Z]/);
	my @line=split /\t/,$_;
	my $keep;
	if($#line>0 && $_=~/[A-z]/){
		$keep=($#line-1)/2;
		print OUT join("\t",@line[0...$keep]),"\n";
		print OUT2 join("\t",@line[$keep+1...$#line]),"\n";
	}else{
		print OUT abs($line[0])+abs($line[1]),"\n";
		print OUT2 abs($line[1]),"\n";
	}
}
`rm $tmp`;
