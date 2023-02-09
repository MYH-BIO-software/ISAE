#!/usr/bin/python3
from optparse import OptionParser
usage = "usage: %prog [options] -s SAE_bed_file -t TAE_bed_file -m Link_gene_module -g  Gene_tss_bed_file  -G Gff_annotation_file  -e Gene_expression_file -o  Output_directory  "
parser = OptionParser(usage = usage)
parser.add_option("-s","--SAE", dest="super_peak_file", default=None,help = "Bed file of SAE.")
parser.add_option("-t","--TAE", dest="typical_peak_file", default=None,help = "Bed file of TAE.")
parser.add_option("-g","--gene_tss", dest="gene_tss_file", default=None,help = "Bed file of gene TSS site.")
parser.add_option("-G","--gff", dest="gff_file", default=None,help = "Gff annotation file of gene.")
parser.add_option("-e","--exp", dest="gene_exp_file", default=None,help = "Gene expression file.")
parser.add_option("-o","--output", dest="out_path", default=None,help = "Path to output directory.")
parser.add_option("-k","--hk", dest="houseKeeping_gene_ID_file", default=None,help = "House keeping gene ID (Optional parameter). If provide this file, these housekeeping genes will be removed in the gene list associated with SAE and TAE.")
parser.add_option("-m","--moudle",dest="link_gene_module", default=None,help = 'The way to link gene with ATAC peak according to the type of SAE and TAE. Could choose "distal" or "proximal"')
(options,args) = parser.parse_args()

if not options.super_peak_file or not options.typical_peak_file or not options.gene_tss_file or not options.gene_exp_file or not options.gff_file  or not options.out_path or not options.link_gene_module  :
    print("Missing parameter, please read the parameter description")
    parser.print_help()
    exit()

if not options.houseKeeping_gene_ID_file :
    delete_HK_gene = "off"
else:
    delete_HK_gene="on"

import os
import subprocess
if options.link_gene_module == "proximal" :
    if os.path.exists(options.out_path) ==False :
        subprocess.call("mkdir "+ options.out_path,shell=True)
    subprocess.call("sort "+ options.gene_exp_file +" > " +options.out_path+"/gene_exp_file_sort" ,shell=True)
    subprocess.call("bedtools intersect -a " +options.super_peak_file+" -b "+ options.gene_tss_file+" -wa -wb > "+ options.out_path + "/super_peak_link_gene_by_tss",shell=True)
    a=subprocess.call("awk '{print $9}' " + options.out_path + "/super_peak_link_gene_by_tss |sort -u > "+ options.out_path+"/gene_ID_tss_link_super_peak_exp" ,shell=True)
    subprocess.call("join "+options.out_path+"/gene_ID_tss_link_super_peak_exp" +" "+ options.out_path+"/gene_exp_file_sort" +" > "+ options.out_path +"/gene_tss_link_super_peak_exp" ,shell=True)
    subprocess.call("bedtools intersect -a " +options.typical_peak_file+" -b "+ options.gene_tss_file+" -wa -wb > "+ options.out_path + "/typical_peak_link_gene_by_tss",shell=True)
    a=subprocess.call("awk '{print $9}' " + options.out_path + "/typical_peak_link_gene_by_tss |sort -u > "+ options.out_path+"/gene_ID_tss_link_typical_peak_exp" ,shell=True)
    subprocess.call("join "+options.out_path+"/gene_ID_tss_link_typical_peak_exp" +" "+ options.out_path+"/gene_exp_file_sort" +" > "+ options.out_path +"/gene_tss_link_typical_peak_exp" ,shell=True)
    subprocess.call(" super_typical_gene_exp_box_plot.R "+options.out_path+"/gene_tss_link_super_peak_exp "+options.out_path+"/gene_tss_link_typical_peak_exp "+options.out_path+"/gene_tss_link_peak_exp_boxplot.tiff",shell=True)
    if delete_HK_gene=="on" :
        subprocess.call("sort "+ options.houseKeeping_gene_ID_file +">" +options.out_path+ "/houseKeeping_gene_ID_file_sort" ,shell=True)
        subprocess.call("join -v1 "+options.out_path+"/gene_tss_link_typical_peak_exp" +" "+ options.out_path+ "/houseKeeping_gene_ID_file_sort > " + options.out_path +"/gene_tss_link_typical_peak_exp_non_houseKeep" ,shell=True)
        subprocess.call("join -v1 "+options.out_path+"/gene_tss_link_super_peak_exp" +" "+ options.out_path+ "/houseKeeping_gene_ID_file_sort > " + options.out_path +"/gene_tss_link_super_peak_exp_non_houseKeep" ,shell=True)
        subprocess.call(" super_typical_gene_exp_box_plot.R "+options.out_path+"/gene_tss_link_super_peak_exp_non_houseKeep "+options.out_path+"/gene_tss_link_typical_peak_exp_non_houseKeep "+options.out_path+"/gene_tss_link_peak_exp_non_houseKeep_boxplot.tiff",shell=True)
elif options.link_gene_module == "distal" :
    if os.path.exists(options.out_path) ==False :
        subprocess.call("mkdir "+ options.out_path,shell=True)
    out_path_distal=options.out_path+"/distal_link_gene"
    if os.path.exists(out_path_distal) ==False :
        subprocess.call("mkdir "+ out_path_distal,shell=True)
    subprocess.call("sort "+ options.gene_exp_file +" > " +options.out_path+"/gene_exp_file_sort" ,shell=True)
    subprocess.call("bedtools intersect -a " +options.super_peak_file+" -b "+ options.gene_tss_file+" -wa -v > "+ out_path_distal + "/super_peak_non_link_gene_by_tss.bed",shell=True)
    subprocess.call("bedtools intersect -a " +options.typical_peak_file+" -b "+ options.gene_tss_file+" -wa -v > "+ out_path_distal + "/typical_peak_non_link_gene_by_tss.bed",shell=True)
    subprocess.call(" chip-seq_annotation.pl -gff " +options.gff_file +" -gffpattern ID -in "+ options.out_path+"  -postfix /*link_gene_by_tss.bed -o "+ options.out_path + " -dis 500000 -id gene --loc both",shell=True)
    subprocess.call(" overlap.pl -in1 " +options.gene_exp_file+ " -in2 " + out_path_distal+ "/super_peak_non_link_gene_by_tss.bed.annoation_500000 -t s -o " + out_path_distal+ "/super_peak_non_link_gene_by_tss.bed.annoation_500000-vs-tmp-gene_name"+ " -s1 '\t' -s2 '\t' -c1 0 -c2 6",shell=True)
    subprocess.call(" overlap.pl -in1 " +options.gene_exp_file+ " -in2 " + out_path_distal+ "/typical_peak_non_link_gene_by_tss.bed.annoation_500000 -t s -o " + out_path_distal+ "/typical_peak_non_link_gene_by_tss.bed.annoation_500000-vs-tmp-gene_name"+ " -s1 '\t' -s2 '\t' -c1 0 -c2 6",shell=True)
    subprocess.call(" pick_chip-seq_annotation_nearest.pl -in "+ options.out_path+" -postfix /*500000-vs-tmp-gene_name -o "+options.out_path+" -loc_index -7 -dis_index -6 -key '\-11,\-10,\-9'",shell=True)
    subprocess.call("awk '{print $1}' "+out_path_distal+"/typical_peak_non_link_gene_by_tss.bed.annoation_500000-vs-tmp-gene_name.gene_nearest  |sort > "+options.out_path+"/gene_ID_distal_link_typical_peak_tmp",shell=True)
    subprocess.call("awk '{print $1}' "+out_path_distal+"/super_peak_non_link_gene_by_tss.bed.annoation_500000-vs-tmp-gene_name.gene_nearest  |sort -u > "+options.out_path+"/gene_ID_distal_link_super_peak",shell=True)
    subprocess.call("join -v1 "+options.out_path+"/gene_ID_distal_link_typical_peak_tmp " +options.out_path+"/gene_ID_distal_link_super_peak |sort -u > "+options.out_path+"/gene_ID_distal_link_typical_peak",shell=True)
    subprocess.call("join "+options.out_path+"/gene_ID_distal_link_typical_peak" +" "+ options.out_path+"/gene_exp_file_sort" +" > "+ options.out_path +"/gene_ID_distal_link_typical_peak_exp",shell=True)
    subprocess.call("join "+options.out_path+"/gene_ID_distal_link_super_peak" +" "+ options.out_path+"/gene_exp_file_sort" +" > "+ options.out_path +"/gene_ID_distal_link_super_peak_exp",shell=True)
    subprocess.call(" super_typical_gene_exp_box_plot.R "+options.out_path+"/gene_ID_distal_link_super_peak_exp "+options.out_path+"/gene_ID_distal_link_typical_peak_exp "+options.out_path+"/gene_distal_link_peak_exp_boxplot.tiff",shell=True)
    if delete_HK_gene=="on" :
        subprocess.call("sort "+ options.houseKeeping_gene_ID_file +">" +options.out_path+ "/houseKeeping_gene_ID_file_sort" ,shell=True)
        subprocess.call("join -v1 "+options.out_path+"/gene_ID_distal_link_super_peak_exp " + options.out_path+ "/houseKeeping_gene_ID_file_sort > "+ options.out_path +"/gene_ID_distal_link_super_peak_exp_non_houseKeep" ,shell=True)
        subprocess.call("join -v1 "+options.out_path+"/gene_ID_distal_link_typical_peak_exp " + options.out_path+ "/houseKeeping_gene_ID_file_sort > " + options.out_path +"/gene_ID_distal_link_typical_peak_exp_non_houseKeep" ,shell=True)
        subprocess.call(" super_typical_gene_exp_box_plot.R "+options.out_path+"/gene_ID_distal_link_super_peak_exp_non_houseKeep "+options.out_path+"/gene_ID_distal_link_typical_peak_exp_non_houseKeep "+options.out_path+"/gene_distal_link_peak_exp_non_houseKeep_boxplot.tiff",shell=True)

