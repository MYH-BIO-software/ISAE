#!/usr/bin/perl -w
my $begin_time=time;
use strict;
use Data::Dumper;
use Getopt::Long;
my ($input,$output,$postfix,$loc,$loc_index,$dis_index,$keep,$keepIndex,$key);
GetOptions(
			"help|?" =>\&USAGE,
			"in:s"=>\$input,
			"o:s"=>\$output,
			"postfix:s"=>\$postfix,
			"loc:s"=>\$loc,
			"loc_index:i"=>\$loc_index,
			"dis_index:i"=>\$dis_index,
			"key:s"=>\$key,
			"keep:s"=>\$keep,
			"keepIndex:i"=>\$keepIndex
			) or &USAGE;
&USAGE unless ($input && $output && $postfix && $loc_index && $dis_index) ;
sub USAGE{
	my $usage=<<"USAGE";
	Description:
	This program can annotation ChIP-seq peak to genes which from gff file.
	Usage:
	-in				input file path					must be given, must be MACS peak file or bed file.
	-postfix				postfix of input files					must be given.
	-o				output file path					must be given. 
	-loc				option(TSS|TES|both)					default is both.
	-loc_index				column index of the loc from the chip-seq annotaion file					must be given.
	-dis_index				column index of the dis from the chip-seq annotaion file					must be given.
	-key				option(default:0,1,2)					The columns which were used as the key of each unit.
	-keep				option(undef)					string which you want to be kept.
	-keepIndex				option(undef)					Which column contain the -keep string you want to be kept.
USAGE
	print $usage;
	exit;
}
####record commend to sh.
$loc="both" unless ($loc);
my @key;
@key=(0,1,2)unless ($key);
my @perl=split/\//,$0;
open SH, ">".$output."/".$perl[-1].".sh" or die $output."$perl[-1].sh\t", $!;
print SH join(" ","perl",$0,"-in",$input,"-postfix",$postfix,"-o",$output,"-loc",$loc,"-loc_index",$loc_index,"-dis_index",$dis_index,"-key",$key)," ";
print SH join(" ","-keep",$keep,"-keepIndex",$keepIndex)if($keep);
print SH "\n";
#####
if($key){
	$key=~s/\\//ig;
	@key=split /,/,$key;
}
print "load file...\nkey is $key \n";
chdir($input);
my @file=glob("*$postfix");
my $add=$loc;
$add="gene"if($loc=~/both/);
for my $each_file(@file){
	my %annotation;
	my @each_file=split /\//,$each_file;
	print $each_file,"\n";
	open OUT,">",$output."/$each_file".".$add\_nearest" or die $output."/$each_file","\t",$!;
	open PEAk, "$input/$each_file" or die $each_file,"\t",$!,"\n";
	my $count=0;
	while(<PEAk>){
		chomp;
		next if(/^@|#/);
		s/\r//ig;
		$count++;
		my @line=split/\t/,$_;
		if($keep){
			next unless ($line[$keepIndex]=~/$keep/);
		}
		my $key2=join("\t",@line[@key]);
		print "The key is: ",$key2,"\n" if($count<=2);
		unless ($loc=~/both/){
			 next unless ($line[$loc_index]=~/$loc/ or $line[$loc_index]=~/internal/);
		}
		push @{$annotation{$key2}},[@line];
		#print $line[$dis_index],"\n";
	}
	my $c=keys %annotation;
	print "in total $c key\n";
	for (keys %annotation){
		my (%dis,$con);
		for (@{$annotation{$_}}){
			my @line=@{$_};
			if($line[$loc_index]=~/internal/){
				print OUT join ("\t",@line),"\n";
				$con=1;
			}else{
				push @{$dis{abs($line[$dis_index])}},join ("\t",@line);
			}
		}
		next if($con);
		my @keys =sort {$a<=>$b} keys %dis;
		for (@{$dis{$keys[0]}}){
			print OUT $_,"\n";
		}
	}
}
