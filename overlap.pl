#!/usr/bin/perl -w
my $begin_time=time;
use strict;
use Data::Dumper;
use Getopt::Long;
my ($input1,$input2,$output,$type,$column1,$column2,$split1,$split2,$split1t,$split2t);
GetOptions(
			"help|?" =>\&USAGE,
			"in1:s"=>\$input1,
			"in2:s"=>\$input2,
			"t:s"=>\$type,
			"s1:s"=>\$split1,
			"s2:s"=>\$split2,
			"st1:s"=>\$split1t,
			"st2:s"=>\$split2t,
			"c1:i"=>\$column1,
			"c2:i"=>\$column2,
			"o:s"=>\$output
				) or &USAGE;
sub USAGE{
	my $usage=<<"USAGE";
	Description:
	This program can find out and statistics overlap or different of first column between two files.
	Usage:
	-in1				input file						input file must be given.
	-in2				input file						input file must be given.
	-t				option(s or d or a)						If choose -t s same lines between 2 files will be output. -t d means output different. -t a means union set.
	-s1				option						If choose this option, input file -in1 will be split as -s1.
	-s2				option						If choose this option, input file -in2 will be split as -s2.
	-c1				option						If  -s was choose -c1 must be given. the column of -in1 that need to be delete redundancy.
	-c2				option						If  -s was choose -c2 must be given. the column of -in2 that need to be delete redundancy.
	-st1				option						split the -in1 line times default: all.
	-st2				option						split the -in2 line times default: all.
	-o				output file						output file must be given.
USAGE
	print $usage;
	exit;
}
&USAGE unless ($input1 && $input2 && $output) ;
$column1=0 if (!$column1 && $split1);
$column2=0 if (!$column2 && $split2);
## get the shell to keep the parameters. 
my @perl=split /\/|\\/,$0; my $sh_output=$output.".sh";
print $perl[-1],"\n";
open SH, ">".$sh_output or die $sh_output,"\t$!";
print SH join(" ","perl",$0,"-in1",$input1,"-in2",$input2,"-t",$type,"-o",$output),"\t";
print SH join(" ","-s1",$split1,"-s2",$split2,"-c1",$column1,"-c2",$column2)if($split2);
print SH "\n";
close SH;
##
my %input1;
open IN, $input1 or die $!;
print "input1: ",$input1,"\n";
my $con=0;
while(<IN>){
	chomp;
	s/\r//ig;
	$con++;
	if($split1){
		my @line;
		if($split1t){
			@line=split (/$split1/,$_,$split1t);
		}else{
			@line=split /$split1/,$_;
		}
		#pop @line unless ($line[-1]);
		pop @line if ($line[-1]=~/^$/);
		push @{$input1{$line[$column1]}},$_;
		#print "aaa$line[$column1]bbb","\n";
		print "xxxxinput1:key $line[$column1] xxxx","\n" if($con<=5);
	}else{
		push @{$input1{$_}},$_;
		print "xxxxinput1:key $_ xxxx","\n" if($con<=5);
		}
}
$con=0;
close IN;
open IN, $input2 or die $!;
print "input2: ",$input2,"\n";
my (%input2,$key);
open OUT, ">".$output or die $!;
print OUT "different line of $input2 to $input1\n"if($type=~/^d$/);
my ($length1,$length2);
while(<IN>){
	chomp;
	$con++;
	s/\r//ig;
	#if($_=~/\t$/)
	if($split2){
		my @line;
		if($split2t){
			@line=split (/$split2/,$_,$split2t);
		}else{
			@line=split /$split2/,$_;
		}
		#print $line[0],"\n" if ($line[-1]=~/^$/);
		pop @line if ($line[-1]=~/^$/);
		$length2=$#line;
		if ($line[$column2]){
			$key =$line[$column2];
		}else{
			print "no $column2 at row $con: ",$_,"\n";
			next
		}
		#print "ccc$line[$column2]ddd","\n";
	}else{
		$key =$_;
	}
	print "xxxxinput2:key $key xxxx","\n" if($con<=5);
	if(!$key && $type=~/^s$/){
		print "here ",$_," omit\n";
		print OUT  $_,"\n";
		next;
	}
	if($type=~/^s$/){
		#print OUT join("\t",@{$input1{$key}},$_),"\n"if(exists $input1{$key});
		&out($input1{$key},\*OUT,$_)if(exists $input1{$key});
	}elsif($type=~/^d$/){
		print OUT $_,"\n" unless(exists $input1{$key});
		$input2{$key}=$_;
	}else{
		$input2{$key}=$_;
		if(exists $input1{$key}){
			#print OUT join("\t",@{$input1{$key}},$_),"\n";
			&out($input1{$key},\*OUT,$_);
			$length1=$#{$input1{$key}};
		}else{
			#$input2{$key}=$_;
			print OUT join("\t",$input1,("not-found")x$length1,$_),"\n";
			#print OUT join("\t",@{$input1{$key}},$_),"\n";
		}
	}
}
#my $input2_t=keys %input1;
#print $input2_t,"\n";
unless($type=~/^s$/){
	print OUT "different line of $input1 to $input2\n" if ($type=~/^d$/);
	for(keys %input1){
		 unless(exists $input2{$_}){
			print OUT join("\t",@{$input1{$_}}),"\n" if ($type=~/^d$/);
			print OUT join("\t",@{$input1{$_}},$input2,("not-found")x$length1),"\n" if ($type=~/^a$/);
			#print join("\t",@{$input1{$_}},$input2,(("\t"."not-found")x$length2)),"\n";
		}
	}
}
print "done\n";
#$a=0;
#print "testing", $a,"\n"if($a=~/^$/);
close IN;
close OUT;

sub out{
	my ($loop,$hanedle,$add)=@_;
	for (@{$loop}){
		if($add){
			print $hanedle join("\t",$_,$add),"\n";
		}else{
			print $hanedle join("\t",$_),"\n";
		}
	}
}