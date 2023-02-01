##!/home/miniconda3/envs/software/bin/python
import getopt, sys

opts, args = getopt.getopt(sys.argv[1:], "g:d:o:h",["pt=","pc=","bt=","bc=","st=","sc=","nt=","nc=","ft=","fc=","fp="])

def usage():
	print( ' python3 idnetify_SOC.py -pt  Treatment ATAC broad peak bed file -pc bam_file  Control ATAC broad bed file -bt  Treatment bam file  -bc  Control bam file  -st  Treatment SAE  bed file  -sc  Control SAE  bed file -nt Treatment sample name -nc Control sample name -o output path \n\n -pt ATAC broad peak bed file of treatment group \n -pc ATAC broad peak bed file of control group \n -bt ATAC bam file of treatment group \n -bc ATAC bam file of control group \n -st SAE  bed file of treatment group -st SAE  bed file of control group \n -nt sample name of treatment group \n -nc sample name of control group \n -o Path to output folder \n-h get help info')

for a,o in opts:
	if a in ('--pt') :
		treatment_peak = o
	elif a in ('--pc') :
		control_peak = o
	elif a in ('--bt') :
		treatment_bam = o
	elif a in ('--bc') :
		control_bam = o
	elif a in ('--st') :
		treatment_SAE = o
	elif a in ('--sc') :
		control_SAE = o
	elif a in ('--nt') :
		treatment_name = o
	elif a in ('--nc') :
		control_name = o
	elif a in ('--fc') :
		control_footprint = o
	elif a in ('--ft') :
		treatment_footprint = o
	elif a in ('--fp') :
		footprint_pvalue = o
	elif a in ('-g') :
		genome = o
	elif a in ('-d') :
		diff_analysis = o
	elif a in ('-o') :
		output_path = o
	elif a in ('-h'):
		usage()
		sys.exit()

import subprocess
import datetime
import os.path
log_file1=open(output_path+"/SAE_differential_analysis.log","a")

##############
treatment_bam_list=treatment_bam.split(",")
control_bam_list=control_bam.split(",")
def bam_makeTag (bam_file_list,out_path,sample_name,sample_list,tag_name) :
    s=0
    for i in bam_file_list:
        s=s+1
#        subprocess.call(' makeTagDirectory '+out_path+"/"+sample_name+"_"+str(s)+" "+i,shell=True)
        tag_name.append(sample_name+"_"+str(s))
        sample_list.append(sample_name)
    return
treatment_sample_list=[]
treatment_tag_list=[]
bam_makeTag(treatment_bam_list,output_path,treatment_name,treatment_sample_list,treatment_tag_list)
control_sample_list=[]
control_tag_list=[]
bam_makeTag(control_bam_list,output_path,control_name,control_sample_list,control_tag_list)
treatment_new_peak=output_path+"/"+treatment_name+"_ATAC_peak.bed"
control_new_peak=output_path+"/"+control_name+"_ATAC_peak.bed"
merge_peak=output_path+"/"+control_name+"_"+treatment_name+"_merge_ATAC_peak.bed"
subprocess.call(" awk '{ print $1 "+'"\t" $2 '+'"\t" $3 "\t" "'+treatment_name+'_"'+" FNR }' "+treatment_peak+" > "+treatment_new_peak,shell=True)
subprocess.call(" awk '{ print $1 "+'"\t" $2 '+'"\t" $3 "\t" "'+control_name+'_"'+" FNR }' "+control_peak+">"+control_new_peak,shell=True)
subprocess.call("cat "+ control_new_peak+ " " + treatment_new_peak+ "| bedtools sort -i - |bedtools merge -i - -c 4 -o distinct > "+merge_peak,shell=True)
Peaks_reads_count=output_path+"/"+treatment_name+"_vs_"+control_name+"_peaks_reads_count"
subprocess.call(' annotatePeaks.pl '+merge_peak+" "+genome+" -size given -raw -d "+ " ".join(control_tag_list)+" "+" ".join(treatment_tag_list) +" 1> " + Peaks_reads_count + " 2> annotatePeaks.log",shell=True)
prefix=treatment_name+"_vs_"+control_name
diff_output=output_path+"/"+prefix+"_diffOutput.txt"
subprocess.call(' getDiffExpression.pl '+Peaks_reads_count+" "+" ".join(control_sample_list)+" "+" ".join(treatment_sample_list)+" -"+diff_analysis+" -export "+prefix+" 1> "+diff_output+" 2> "+output_path+"/getDiffExpression.log",shell=True)
signal_up_peak=output_path+"/"+prefix+"_peak_signal_up_FC_1.5_FDR_0.01.txt"
signal_down_peak=output_path+"/"+prefix+"_peak_signal_down_FC_1.5_FDR_0.01.txt"
subprocess.call(' awk -F "\t" '+ "  ' $24 >= 0.58 && $26<0.01  {print $1 "+'"\t" $2 '+'"\t" $3 "\t" $4'+'"\t" $24 "\t" $26'+"  }' "+diff_output+" | sort -k 1 > "+signal_up_peak,shell=True)
subprocess.call(' awk -F "\t" '+ "  ' $24 <= -0.58 && $26<0.01  {print $1 "+'"\t" $2 '+'"\t" $3 "\t" $4'+'"\t" $24 "\t" $26'+"  }' "+diff_output+" | sort -k 1 > "+signal_down_peak,shell=True)
###############
treatment_footprint_pfilter=output_path+"/"+treatment_name+"_footprint_pfilter.bed"
control_footprint_pfilter=output_path+"/"+control_name+"_footprint_pfilter.bed"
def footprint_Num_stat (footprint_file,peak_file,out_path,sample_name,pvalue) :
    footprint_pfilter=out_path+"/"+sample_name+"_footprint_pfilter.bed"
    subprocess.call(" awk ' $5 <  " + pvalue + " ' "+footprint_file+" > "+footprint_pfilter,shell=True)
    peak_number_1=out_path+"/"+sample_name+"_more_than_0_footprint_number"
    peak_number_2=out_path+"/"+sample_name+"_0_footprint_number"
    subprocess.call(" bedtools intersect -a   " + peak_file + " -b "+footprint_pfilter+" -wa -wb | awk '{print $4}' |sort |uniq -c |awk '{print $2 "+' "\t" '+ " $1}' > "+ peak_number_1,shell=True)
    subprocess.call(" bedtools intersect -a   " + peak_file + " -b "+footprint_pfilter+" -wa -v | awk '{print $4  "+' "\t" '+ " 0 }' > "+ peak_number_2,shell=True)
    subprocess.call(" cat " +peak_number_1 + " " + peak_number_2 + " > " + out_path+"/"+sample_name+"_footprint_number",shell=True)
    return
footprint_Num_stat(treatment_footprint,treatment_new_peak,output_path,treatment_name,footprint_pvalue)
footprint_Num_stat(control_footprint,control_new_peak,output_path,control_name,footprint_pvalue)
treatment_footprint_Num=output_path+"/"+treatment_name+"_footprint_number"
control_footprint_Num=output_path+"/"+control_name+"_footprint_number"
Num_change_output=output_path+"/"+prefix+"_footprint_change.txt"
subprocess.call(" identify_merge_broadpeak_footprint_number.py  -t " +treatment_footprint_Num + " -c " + control_footprint_Num + " --nt " + treatment_name +" --nc "+control_name+ " -i " +merge_peak + " -o " + Num_change_output,shell=True)
TF_up_peak=output_path+"/"+prefix+"_peak_footprint_up.txt"
TF_down_peak=output_path+"/"+prefix+"_peak_footprint_down.txt"
subprocess.call(" awk ' $4 > 0 ' " +Num_change_output + " | sort -k 1 > " + TF_up_peak,shell=True)
subprocess.call(" awk ' $4 < 0 ' " +Num_change_output + " | sort -k 1 > " + TF_down_peak,shell=True)
peak_signal_TF_up=output_path+"/"+prefix+"_peak_signal_footprint_up.bed"
peak_signal_TF_down=output_path+"/"+prefix+"_peak_signal_footprint_down.bed"

subprocess.call(" signal_footprint_change.py  -s  " +signal_up_peak + " -f " +TF_up_peak + " -o " + peak_signal_TF_up,shell=True)
subprocess.call(" signal_footprint_change.py  -s  " +signal_down_peak+ " -f " + TF_down_peak + " -o " + peak_signal_TF_down,shell=True)
up_SAE=output_path+"/"+prefix+"_up_SAE.bed"
subprocess.call(" bedtools intersect -a   " + peak_signal_TF_up + " -b "+treatment_SAE+" -wa -u > "+ up_SAE,shell=True)
down_SAE=output_path+"/"+prefix+"_down_SAE.bed"
subprocess.call(" bedtools intersect -a   " + peak_signal_TF_down + " -b "+control_SAE+" -wa -u > "+ down_SAE,shell=True)
subprocess.call(" rm   " + output_path+"/"+"*_more_than_0_footprint_number " + output_path+"/"+"*_0_footprint_number " ,shell=True)
