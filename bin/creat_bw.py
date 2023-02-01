#!/usr/bin/python3
import getopt, sys

opts, args = getopt.getopt(sys.argv[1:], "b:g:n:Ah",["region=","outdir="])

def usage():
	print( ' --region ATAC_broad_peak   -b ATAC_bam_file\n -A ATAC moudle -g chrom sizes file \n -n Experiment name, which will be used to generate output file names  \n --outdir OUTDIR\n -h get help info')

ATAC_module =' ' 
import os
default_path=os.getcwd()
output_path = default_path
for a,o in opts:
	if a in ('--region'):
		broad_open_chromatin_file = o
	if a in ('-g'):
		chrom_size_file = o
	if a in ('-n'):
		output_name = o
	if a in ('-b'):
		input_bam_file = o
	if a in ('--outdir'):
		output_path = o
	if a in ('-A'):
		ATAC_module ="-A" 
	if a in ('-h'):
		usage()
		sys.exit()
import subprocess
creat_dnase_wig_command="dnase_wig_tracks.py "+ broad_open_chromatin_file +" " +input_bam_file+" "+ output_path+"/"+output_name+"_danse_fr.wig "+" "+ output_path+"/"+output_name+"_danse_rev.wig "+ATAC_module
subprocess.call(creat_dnase_wig_command,shell=True)
merge_fr_rev_wig_command="perl encode_ATAC-seq_wig_merge.pl -fr " +output_path+"/"+output_name+"_danse_fr.wig -rev " +  output_path+"/"+output_name+"_danse_rev.wig -out " +output_path+"/"+output_name+"_danse_merge_fr_rev.wig"
subprocess.call(merge_fr_rev_wig_command,shell=True)
subprocess.call("wigToBigWig -clip "+output_path+"/"+output_name+"_danse_merge_fr_rev.wig " + chrom_size_file +" "+output_path+"/"+output_name+"_danse_merge_fr_rev.bw",shell=True)
