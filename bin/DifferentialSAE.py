#!/usr/bin/python3
from optparse import OptionParser
usage = "usage: %prog [options]  --pt  Treatment_ATAC_broad_peak_bed_file --pc   Control_ATAC_broad_bed_file --bt  Treatment_bam_file  --bc  Control_bam_file  --st  Treatment_SAE_bed_file  --sc  Control_SAE_bed_file --nt Treatment_sample_name --nc Control_sample_name --ft Treatment_footprint_peak_bed_file  --fc Control_footprint_bed_file --fs Footprint_Score_threshold -g Genome -d Diff_analysis_method -o output_path   "
parser = OptionParser(usage = usage)
parser.add_option("--pt", dest="treatment_peak", default=None,help = "ATAC broad peak bed file of treatment group.")
parser.add_option("--pc", dest="control_peak", default=None,help = "ATAC broad peak bed file of control group.")
parser.add_option("--bt", dest="treatment_bam", default=None,help = "ATAC bam file of treatment group.")
parser.add_option("--bc", dest="control_bam", default=None,help = "ATAC bam file of control group.")
parser.add_option("--st", dest="treatment_SAE", default=None,help = "SAE bed file of treatment group.")
parser.add_option("--sc", dest="control_SAE", default=None,help = "TAE bed file of control group.")
parser.add_option("--nt", dest="treatment_name", default=None,help = "Experiment name of treatment group.")
parser.add_option("--nc", dest="control_name", default=None,help = "Experiment name of control group.")
parser.add_option("--ft", dest="treatment_footprint", default=None,help = "Footprint bed file of treatment group.")
parser.add_option("--fc", dest="control_footprint", default=None,help = "Footprint bed file of control group.")
parser.add_option("--fs", dest="footprint_score", default=None,help = "Score threshold of footprint of treatment and control group. The footprint with score > set threshold would be used to caclulate footprint number in ATAC broad peak. Recommending minimum footprint score threshold for use in SAE identification of treatment and control group.")
parser.add_option("-g", dest="genome", default=None,help = "Available Genomes:[mm10,hg38,oviAri4,susScr11,bosTau9,galGal5],If no genome is available, specify 'none'.")
parser.add_option("-d", dest="diff_analysis", default="edgeR",help = "Differential Expression program selection:[DESeq2,DESeq,edgeR,limma], Default:edgeR.")
parser.add_option("-o", dest="output_path", default=None,help = "Path to output directory.")
(options,args) = parser.parse_args()

if not options.treatment_peak or not options.control_peak or not options.treatment_bam or not options.control_bam or not options.treatment_SAE  or not  options.control_SAE   or not   options.treatment_name  or not  options.control_name   or not   options.treatment_footprint  or not  options.control_footprint    or not options.footprint_score    or not   options.genome or not  options.output_path :
    print("Missing parameter, please read the parameter description")
    parser.print_help()
    exit()

import subprocess
import datetime
import os.path
log_file1=open(options.output_path+"/SAE_differential_analysis.log","a")

##############
treatment_bam_list=options.treatment_bam.split(",")
control_bam_list=options.control_bam.split(",")
def bam_makeTag (bam_file_list,out_path,sample_name,sample_list,tag_name) :
    s=0
    for i in bam_file_list:
        s=s+1
        subprocess.call(' makeTagDirectory '+out_path+"/"+sample_name+"_"+str(s)+" "+i,shell=True)
        tag_name.append(sample_name+"_"+str(s))
        sample_list.append(sample_name)
    return
treatment_sample_list=[]
treatment_tag_list=[]
bam_makeTag(treatment_bam_list,options.output_path,options.treatment_name,treatment_sample_list,treatment_tag_list)
control_sample_list=[]
control_tag_list=[]
bam_makeTag(control_bam_list,options.output_path,options.control_name,control_sample_list,control_tag_list)
treatment_new_peak=options.output_path+"/"+options.treatment_name+"_ATAC_peak.bed"
control_new_peak=options.output_path+"/"+options.control_name+"_ATAC_peak.bed"
merge_peak=options.output_path+"/"+options.control_name+"_"+options.treatment_name+"_merge_ATAC_peak.bed"
subprocess.call(" awk '{ print $1 "+'"\t" $2 '+'"\t" $3 "\t" "'+options.treatment_name+'_"'+" FNR }' "+options.treatment_peak+" > "+treatment_new_peak,shell=True)
subprocess.call(" awk '{ print $1 "+'"\t" $2 '+'"\t" $3 "\t" "'+options.control_name+'_"'+" FNR }' "+options.control_peak+">"+control_new_peak,shell=True)
subprocess.call("cat "+ control_new_peak+ " " + treatment_new_peak+ "| bedtools sort -i - |bedtools merge -i - -c 4 -o distinct > "+merge_peak,shell=True)
Peaks_reads_count=options.output_path+"/"+options.treatment_name+"_vs_"+options.control_name+"_peaks_reads_count"
subprocess.call(' annotatePeaks.pl '+merge_peak+" "+options.genome+" -size given -raw -d "+ " ".join(control_tag_list)+" "+" ".join(treatment_tag_list) +" 1> " + Peaks_reads_count + " 2> annotatePeaks.log",shell=True)
prefix=options.treatment_name+"_vs_"+options.control_name
diff_output=options.output_path+"/"+prefix+"_diffOutput.txt"
subprocess.call(' getDiffExpression.pl '+Peaks_reads_count+" "+" ".join(control_sample_list)+" "+" ".join(treatment_sample_list)+" -"+options.diff_analysis+" -export "+prefix+" 1> "+diff_output+" 2> "+options.output_path+"/getDiffExpression.log",shell=True)
signal_up_peak=options.output_path+"/"+prefix+"_peak_signal_up_FC_1.5_FDR_0.01.txt"
signal_down_peak=options.output_path+"/"+prefix+"_peak_signal_down_FC_1.5_FDR_0.01.txt"
subprocess.call(' awk -F "\t" '+ "  ' $24 >= 0.58 && $26<0.01  {print $1 "+'"\t" $2 '+'"\t" $3 "\t" $4'+'"\t" $24 "\t" $26'+"  }' "+diff_output+" | sort -k 1 > "+signal_up_peak,shell=True)
subprocess.call(' awk -F "\t" '+ "  ' $24 <= -0.58 && $26<0.01  {print $1 "+'"\t" $2 '+'"\t" $3 "\t" $4'+'"\t" $24 "\t" $26'+"  }' "+diff_output+" | sort -k 1 > "+signal_down_peak,shell=True)
###############
treatment_footprint_pfilter=options.output_path+"/"+options.treatment_name+"_footprint_pfilter.bed"
control_footprint_pfilter=options.output_path+"/"+options.control_name+"_footprint_pfilter.bed"
def footprint_Num_stat (footprint_file,peak_file,out_path,sample_name,score) :
    footprint_sfilter=out_path+"/"+sample_name+"_footprint_sfilter.bed"
    subprocess.call(" awk ' $5 <  " + score + " ' "+footprint_file+" > "+footprint_sfilter,shell=True)
    peak_number_1=out_path+"/"+sample_name+"_more_than_0_footprint_number"
    peak_number_2=out_path+"/"+sample_name+"_0_footprint_number"
    subprocess.call(" bedtools intersect -a   " + peak_file + " -b "+footprint_sfilter+" -wa -wb | awk '{print $4}' |sort |uniq -c |awk '{print $2 "+' "\t" '+ " $1}' > "+ peak_number_1,shell=True)
    subprocess.call(" bedtools intersect -a   " + peak_file + " -b "+footprint_sfilter+" -wa -v | awk '{print $4  "+' "\t" '+ " 0 }' > "+ peak_number_2,shell=True)
    subprocess.call(" cat " +peak_number_1 + " " + peak_number_2 + " > " + out_path+"/"+sample_name+"_footprint_number",shell=True)
    return
footprint_Num_stat(options.treatment_footprint,treatment_new_peak,options.output_path,options.treatment_name,options.footprint_score)
footprint_Num_stat(options.control_footprint,control_new_peak,options.output_path,options.control_name,options.footprint_score)
treatment_footprint_Num=options.output_path+"/"+options.treatment_name+"_footprint_number"
control_footprint_Num=options.output_path+"/"+options.control_name+"_footprint_number"
Num_change_output=options.output_path+"/"+prefix+"_footprint_change.txt"
subprocess.call(" identify_merge_broadpeak_footprint_number.py  -t " +treatment_footprint_Num + " -c " + control_footprint_Num + " --nt " + options.treatment_name +" --nc "+options.control_name+ " -i " +merge_peak + " -o " + Num_change_output,shell=True)
TF_up_peak=options.output_path+"/"+prefix+"_peak_footprint_up.txt"
TF_down_peak=options.output_path+"/"+prefix+"_peak_footprint_down.txt"
subprocess.call(" awk ' $4 > 0 ' " +Num_change_output + " | sort -k 1 > " + TF_up_peak,shell=True)
subprocess.call(" awk ' $4 < 0 ' " +Num_change_output + " | sort -k 1 > " + TF_down_peak,shell=True)
peak_signal_TF_up=options.output_path+"/"+prefix+"_peak_signal_footprint_up.bed"
peak_signal_TF_down=options.output_path+"/"+prefix+"_peak_signal_footprint_down.bed"

subprocess.call(" signal_footprint_change.py  -s  " +signal_up_peak + " -f " +TF_up_peak + " -o " + peak_signal_TF_up,shell=True)
subprocess.call(" signal_footprint_change.py  -s  " +signal_down_peak+ " -f " + TF_down_peak + " -o " + peak_signal_TF_down,shell=True)
up_SAE=options.output_path+"/"+prefix+"_up_SAE.bed"
subprocess.call(" bedtools intersect -a   " + peak_signal_TF_up + " -b "+options.treatment_SAE+" -wa -u > "+ up_SAE,shell=True)
down_SAE=options.output_path+"/"+prefix+"_down_SAE.bed"
subprocess.call(" bedtools intersect -a   " + peak_signal_TF_down + " -b "+options.control_SAE+" -wa -u > "+ down_SAE,shell=True)
subprocess.call(" rm   " + options.output_path+"/"+"*_more_than_0_footprint_number " + options.output_path+"/"+"*_0_footprint_number " ,shell=True)
