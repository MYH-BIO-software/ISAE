#!/home/miniconda3/envs/software/bin/python3
import getopt, sys

opts, args = getopt.getopt(sys.argv[1:], "s:t:g:G:e:o:k:h")

def usage():
	print( ' -s super peak bed file \n -t typical peak file\n -g  gene tss bed file \n -G gff annotation file\n -k houseKeepgene_ID_list -h get help info')

delete_HK_gene="off"
import os
default_path=os.getcwd()
out_path = default_path
for a,o in opts:
	if a in ('-s'):
		super_peak_file = o
	if a in ('-t'):
		typical_peak_file = o
	if a in ('-g'):
		gene_tss_file  = o
	if a in ('-e'):
		gene_exp_file  = o
	if a in ('-o'):
		out_path = o
	if a in ('-G'):
		gff_file = o
	if a in ('-k'):
		houseKeeping_gene_ID_file = o
		delete_HK_gene="on"
	if a in ('-h'):
		usage()
		sys.exit()
import subprocess
subprocess.call("sort "+ gene_exp_file +" > " +out_path+"/gene_exp_file_sort" ,shell=True)
subprocess.call("bedtools intersect -a " +super_peak_file+" -b "+ gene_tss_file+" -wa -wb > "+ out_path + "/super_peak_link_gene_by_tss",shell=True)
a=subprocess.call("awk '{print $9}' " + out_path + "/super_peak_link_gene_by_tss |sort -u > "+ out_path+"/gene_ID_tss_link_super_peak_exp" ,shell=True)
subprocess.call("join "+out_path+"/gene_ID_tss_link_super_peak_exp" +" "+ out_path+"/gene_exp_file_sort" +" > "+ out_path +"/gene_tss_link_super_peak_exp" ,shell=True)
subprocess.call("bedtools intersect -a " +typical_peak_file+" -b "+ gene_tss_file+" -wa -wb > "+ out_path + "/typical_peak_link_gene_by_tss",shell=True)
a=subprocess.call("awk '{print $9}' " + out_path + "/typical_peak_link_gene_by_tss |sort -u > "+ out_path+"/gene_ID_tss_link_typical_peak_exp" ,shell=True)
subprocess.call("join "+out_path+"/gene_ID_tss_link_typical_peak_exp" +" "+ out_path+"/gene_exp_file_sort" +" > "+ out_path +"/gene_tss_link_typical_peak_exp" ,shell=True)
if delete_HK_gene=="on" :
    subprocess.call("sort "+ houseKeeping_gene_ID_file +">" +out_path+ "/houseKeeping_gene_ID_file_sort" ,shell=True)
    subprocess.call("join -v1 "+out_path+"/gene_tss_link_typical_peak_exp" +" "+ out_path+ "/houseKeeping_gene_ID_file_sort > " + out_path +"/gene_tss_link_typical_peak_exp_non_houseKeep" ,shell=True)
    subprocess.call("join -v1 "+out_path+"/gene_tss_link_super_peak_exp" +" "+ out_path+ "/houseKeeping_gene_ID_file_sort > " + out_path +"/gene_tss_link_super_peak_exp_non_houseKeep" ,shell=True)
if  os.path.exists(out_path +"/distal_interaction") ==False :
    subprocess.call("mkdir "+ out_path +"/distal_interaction" ,shell=True)
out_path_distal=out_path +"/distal_interaction" 
subprocess.call("bedtools intersect -a " +super_peak_file+" -b "+ gene_tss_file+" -wa -v > "+ out_path_distal + "/super_peak_non_link_gene_by_tss.bed",shell=True)
subprocess.call("bedtools intersect -a " +typical_peak_file+" -b "+ gene_tss_file+" -wa -v > "+ out_path_distal + "/typical_peak_non_link_gene_by_tss.bed",shell=True)
subprocess.call(" chip-seq_annotation.pl -gff " +gff_file +" -gffpattern ID -in "+ out_path+"  -postfix /*link_gene_by_tss.bed -o "+ out_path + " -dis 500000 -id gene --loc both",shell=True)
subprocess.call(" overlap.pl -in1 " +gene_exp_file+ " -in2 " + out_path_distal+ "/super_peak_non_link_gene_by_tss.bed.annoation_500000 -t s -o " + out_path_distal+ "/super_peak_non_link_gene_by_tss.bed.annoation_500000-vs-tmp-gene_name"+ " -s1 '\t' -s2 '\t' -c1 0 -c2 6",shell=True)
subprocess.call(" overlap.pl -in1 " +gene_exp_file+ " -in2 " + out_path_distal+ "/typical_peak_non_link_gene_by_tss.bed.annoation_500000 -t s -o " + out_path_distal+ "/typical_peak_non_link_gene_by_tss.bed.annoation_500000-vs-tmp-gene_name"+ " -s1 '\t' -s2 '\t' -c1 0 -c2 6",shell=True)
subprocess.call(" pick_chip-seq_annotation_nearest.pl -in "+ out_path+" -postfix /*500000-vs-tmp-gene_name -o "+out_path+" -loc_index -7 -dis_index -6 -key '\-11,\-10,\-9'",shell=True)
subprocess.call("awk '{print $1}' "+out_path_distal+"/typical_peak_non_link_gene_by_tss.bed.annoation_500000-vs-tmp-gene_name.gene_nearest  |sort > "+out_path_distal+"/gene_ID_distal_link_typical_peak_tmp",shell=True)
subprocess.call("awk '{print $1}' "+out_path_distal+"/super_peak_non_link_gene_by_tss.bed.annoation_500000-vs-tmp-gene_name.gene_nearest  |sort -u > "+out_path_distal+"/gene_ID_distal_link_super_peak",shell=True)
subprocess.call("join -v1 "+out_path_distal+"/gene_ID_distal_link_typical_peak_tmp " +out_path_distal+"/gene_ID_distal_link_super_peak |sort -u > "+out_path_distal+"/gene_ID_distal_link_typical_peak",shell=True)
subprocess.call("join "+out_path_distal+"/gene_ID_distal_link_typical_peak" +" "+ out_path+"/gene_exp_file_sort" +" > "+ out_path_distal +"/gene_ID_distal_link_typical_peak_exp",shell=True)
subprocess.call("join "+out_path_distal+"/gene_ID_distal_link_super_peak" +" "+ out_path+"/gene_exp_file_sort" +" > "+ out_path_distal +"/gene_ID_distal_link_super_peak_exp",shell=True)
subprocess.call(" super_typical_gene_exp_box_plot.R "+out_path+"/gene_tss_link_super_peak_exp "+out_path+"/gene_tss_link_typical_peak_exp "+out_path+"/gene_tss_link_peak_exp_boxplot.tiff",shell=True)
subprocess.call(" super_typical_gene_exp_box_plot.R "+out_path_distal+"/gene_ID_distal_link_super_peak_exp "+out_path_distal+"/gene_ID_distal_link_typical_peak_exp "+out_path_distal+"/gene_distal_link_peak_exp_boxplot.tiff",shell=True)
if delete_HK_gene=="on" :
    subprocess.call("join -v1 "+out_path_distal+"/gene_ID_distal_link_super_peak_exp " + out_path+ "/houseKeeping_gene_ID_file_sort > "+ out_path_distal +"/gene_ID_distal_link_super_peak_exp_non_houseKeep" ,shell=True)
    subprocess.call("join -v1 "+out_path_distal+"/gene_ID_distal_link_typical_peak_exp " + out_path+ "/houseKeeping_gene_ID_file_sort > " + out_path_distal +"/gene_ID_distal_link_typical_peak_exp_non_houseKeep" ,shell=True)
    subprocess.call(" super_typical_gene_exp_box_plot.R "+out_path+"/gene_tss_link_super_peak_exp_non_houseKeep "+out_path+"/gene_tss_link_typical_peak_exp_non_houseKeep "+out_path+"/gene_tss_link_peak_exp_non_houseKeep_boxplot.tiff",shell=True)
    subprocess.call(" super_typical_gene_exp_box_plot.R "+out_path_distal+"/gene_ID_distal_link_super_peak_exp_non_houseKeep "+out_path_distal+"/gene_ID_distal_link_typical_peak_exp_non_houseKeep "+out_path_distal+"/gene_distal_link_peak_exp_non_houseKeep_boxplot.tiff",shell=True)

