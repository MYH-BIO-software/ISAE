#!/usr/bin/perl -w
my $begin_time=time;
use strict;
use Data::Dumper;
use Getopt::Long;
my ($gff_file,$input,$output,$postfix,$loc,$intergenic_out,$dis,$submit,$gffpattern,$transcript_ID);
GetOptions(
			"help|?" =>\&USAGE,
			"gff:s"=>\$gff_file,
			"gffpattern:s"=>\$gffpattern,
			"in:s"=>\$input,
			"o:s"=>\$output,
			"postfix:s"=>\$postfix,
			"loc:s"=>\$loc,
			"i"=>\$intergenic_out,
			"dis:s"=>\$dis,
			"submit:s"=>\$submit,
			"id:s"=>\$transcript_ID,
			) or &USAGE;
&USAGE unless ($gff_file && $input && $output && $postfix ) ;
sub USAGE{
	my $usage=<<"USAGE";
	Description:
	This program can annotation ChIP-seq peak to genes which from gff file.
	####the postive number of dis is the peak located at the upstream of setted -loc of gene and the minus is the downstream.
	Usage:
	-gff				GFF file					must be given.
	-gffpattern				option(defaut:gene_name)					The pattern which were used to matched the gff file of -id line in the program 
	-in				input file path					must be given, must be MACS peak file or bed file.
	-postfix				postfix of input files					must be given.
	-o				output file path					must be given. 
	-loc				option(TSS|TES|both)					default is both. 
	-i				option(default:undeaf)					if choose this option intergenic peak will be output.
	-dis				option(default:5000)					Must be digit. The distance between peak and gene will be associated.
	--submit				option(default:undeaf)					If this option was given the --submit column will use as the peak location. If undeaf the cent of peak will use as peak location.
	-id				option(default:exon)					The pattern which were used to matched the gff file which will be use as the transcript information.
USAGE
	print $usage;
	exit;
}
####record commend to sh.
$transcript_ID="exon"unless ($transcript_ID);
$loc ="both"unless ($loc );
my @perl=split/\//,$0;
open SH, ">".$output."/".$perl[-1].".sh" or die $output."$perl[-1].sh\t", $!;
print SH "#the postive number of dis is the peak located at the upstream of expected $loc loc of gene and the minus is the downstream.\n";
print SH join(" ","perl",$0,"-gff",$gff_file,"-gffpattern",$gffpattern,"-in",$input,"-postfix",$postfix,"-o",$output,"-dis",$dis,"-id",$transcript_ID)," ";
print SH join(" ","--submit",$submit)if($submit);
print SH join(" ","--loc",$loc)if($loc);
print SH "\n";
#####
my (%ori_chr,%chr);
my $base=50000;
print "loads gff...\n";
open GFF,$gff_file or die $!;
my (%exon,$key);
while(<GFF>){
	chomp;
	next if($_=~/^$|#/);
	my @line=split /\t/,$_;
	next unless (@line);
	if($line[2]=~/gene/){
		$key=$_;
		push @{$ori_chr{$line[0]}},$_;
		#print $_,"\n";
	}
	if($line[2]=~/$transcript_ID/){###keep exon ingormation.
		#push @{$exon{$key}},[@line];
		push @{$exon{$key}},[@line];
	}
	#print $_ ,"\n"if($_=~/LOC100621233/);
}
close GFF;
print "loads gff done.\nblock gff....\n";
my ($start,$end,@block);
for my $keys(keys %ori_chr){
	#print "$keys\n";
	for my $line (@{$ori_chr{$keys}}){
		my @line=split /\t/,$line;
		$start=int($line[3]/$base);
		$end=int($line[4]/$base);
		for ($start...$end){
			push @{$block[$_]},$line;
			print join("\t",@line),"\t","XLOC_045578 index ",$_,"\n"if($line[-1]=~/XLOC_045578/);
		}
	}
	push @{$chr{$keys}},@block;
	#print Dumper $block[$start],"\n";
	@block=();
}
undef %ori_chr;
#print Dumper @{${$chr{'chr11'}}[323]},"\n";
print "block gff done.\nload peak file...\n";
chdir($input);
my @file=glob("*$postfix");
#print @file,"\n";
for my $each_file(@file){
	my @each_file=split /\//,$each_file;
	print $each_file,"\n";
	open INTER, ">$each_file\_intergenic_peak.bed" or die $! if($intergenic_out);
	open OUT,">",$output."/$each_file".".annoation_$dis" or die $output."/$each_file","\t",$!;
	open PEAk, "$input/$each_file" or die $each_file,"\t",$!,"\n";
	while(<PEAk>){
		chomp;
		next if(/^@/);
		s/\r//ig;
		my @line=split/\t/,$_;
		if ($submit && $line[$submit]=~/[A-z]/){
			print OUT join("\t",@line),"\n";
			print join("\t",@line),"omit\n";
			next 
		}elsif ($line[1]=~/[A-z]/||$line[2]=~/[A-z]/){
			print OUT join("\t",@line),"\n";
			print join("\t",@line),"omit\n";
			next
		}
		my ($peak_start,$peak_end,$peak);
		if ($submit){
			$peak_start=$line[$submit]-$dis;
			$peak_end=$line[$submit]+$dis;
			$peak=$line[$submit];
			$peak_start=1 if($peak_start<0);
		}else{
			$peak=($line[1]+$line[2])/2;
			$peak_start=$peak-$dis;
			$peak_end=$peak+$dis;
			$peak_start=1 if($peak_start<0);
		}
		#print $dis,"\n";
		if(exists $chr{$line[0]}){#需要考虑gff 和 reference seq之间对应不上：eg：chrMT在gff中就没有。
			my $start=int($peak_start/$base);
			$end=int($peak_end/$base);#要考虑长度问题，比如太长将基因包进去，只看start 和 end 时就考虑不到。
			my @loop;
			for($start...$end){
				push @loop,@{$chr{$line[0]}[$_]} if($chr{$line[0]}[$_]);
			}
			unless (@loop){#当reads的起始点没有落在染色体分块中时检查reads结束点是否落在某一染色体分块中
				#$intergenic++;
				print INTER join("\t",@line),"\n"if($intergenic_out);
				next;
			}
			#print join("\t",@line,$peak_start,$peak_end,$start,$end,$base),"\n",Dumper @loop,"\n" if($line[3]=~/3_15802_lociStitched/);
			my (%output);
			for (@loop){
				#next unless(@{$_});
				my @gene_line=split /\t/,$_;
				#print $_,"\n";
				next if($gene_line[3]>$peak_end or $gene_line[4]<$peak_start);# 检查peak是否落在某个基因内部
				######
				#my $gene_name=0;
				my (%gene_name,$split);
				for(@{$exon{$_}}){
					my @exon=@{$_};
					if($_=~/";/){
						$split="\";"
					}else{
						$split=";"
					}
					my @gene_name=split/$split/,$exon[8];
					#print join("\t",@gene_name),"\n" if($line[3]=~/3_15802_lociStitched/);
					for(@gene_name){
						$gene_name{$2}=1 if($_=~/$gffpattern( "|\=)(\w.*)/);
						#print "ok\n" if($_=~/Name/);
						#print $_,$1,"aaaa \n" if($_=~/$gffpattern\=(\w.*)/);;
					}
				}
				#$gene_name=join (",",keys %gene_name) if(%gene_name);
				#unless ( $gene_name){
				unless ( %gene_name){
					my @gene_name=split/$split/,$gene_line[8];
					for(@gene_name){
						#$gene_name=$2 if($_=~/(Name|$gffpattern)\=(\w.*)/);
						$gene_name{$2}=1 if($_=~/(Name|$gffpattern)\=(\w.*)/);
					}
				}
				#print $gene_name,"\n";
				#######
				my ($strand_loc,$dis2,$dis_s,$dis_e);
				if($gene_line[3]<$peak && $peak<$gene_line[4]){
					#print OUT,join("\t",@line,"internal","0",$gene_name,@gene_line[0...2,8]);
					$dis_s=$peak-$gene_line[3];
					$dis_e=$peak-$gene_line[4];
					if($gene_line[6]=~/\Q+\E/){
						if($loc=~/TSS/ && abs($dis_s)<$dis){
							$strand_loc="internal_TSS";
							$dis2=$dis_s;
						}elsif($loc=~/TES/ && abs($dis_e)<$dis){
							$strand_loc="internal_TES";
							$dis2=$dis_e;
						}elsif(abs($dis_e)>abs($dis_s)){
							$strand_loc="internal_TSS";
							$dis2=$dis_s;
						}else{
							$strand_loc="internal_TES";
							$dis2=$dis_e;
						}
					}else{
						if($loc=~/TSS/ && abs($dis_s)<$dis){
							$strand_loc="internal_TSS";
							$dis2=-$dis_e;
						}elsif($loc=~/TES/ && abs($dis_e)<$dis){
							$strand_loc="internal_TES";
							$dis2=-$dis_s;
						}elsif(abs($dis_e)>abs($dis_s)){
							$strand_loc="internal_TES";
							$dis2=-$dis_s;
						}else{
							$strand_loc="internal_TSS";
							$dis2=-$dis_e;
						}
					}
				}elsif(0<=$gene_line[3]-$peak){
					if($gene_line[6]=~/\Q+\E/){
						$strand_loc="TSS";
						$dis2=$peak-$gene_line[3];
					}else{
						$strand_loc="TES";
						$dis2=$gene_line[3]-$peak;
					}
				}else{
					if($gene_line[6]=~/\Q+\E/){
						$strand_loc="TES";
						$dis2=$peak-$gene_line[4];
					}else{
						$strand_loc="TSS";
						$dis2=$gene_line[4]-$peak
						}
				}
				#$gene_name="NA"unless ($gene_name);
				$gene_name{"NA"}=1 unless (%gene_name);
				#if($loc=~/TSS/ && $strand_loc=~/TSS/ && abs($dis_s)<$dis){$strand_loc="internal_TSS"; $dis2=abs($dis_s)}
				if($loc=~/both/){
					#print OUT join("\t",@line,$strand_loc,$dis2,$gene_name,@gene_line[0,3,4,8]),"\n";
					&out(\%output,\%gene_name,@line,$strand_loc,$dis2,"NA",@gene_line[0,3,4,8]);
					
				}elsif ($strand_loc=~/$loc/){
					#print OUT join("\t",@line,$strand_loc,$dis2,$gene_name,@gene_line[0,3,4,8]),"\n";
					&out(\%output,\%gene_name,@line,$strand_loc,$dis2,"NA",@gene_line[0,3,4,8]);
				}
			}
			if(%output){
				for (keys %output){
					print OUT $_,"\n";
				}
			}
		}
	}
}

sub out{
	my $output=shift @_;
	my $loop=shift @_;
	my @print=@_;
	for (keys %{$loop}){
		$print[-5]=$_;
		${$output}{ join("\t",@print)}=1;
	}
}