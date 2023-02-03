#!/usr/bin/python3
#将ATAC peak区分为proximal和distal两类分别鉴定，并可以选择除去落在houseKeeping基因TSS region的ATAC peak
import getopt, sys

opts, args = getopt.getopt(sys.argv[1:], "t:f:b:p:g:n:Ah",["outdir=","g_file=","Pvalue=","macs2_call_peak=","wellington_footprints=","Broad_peak=","Footprint_file=","Dnase_wig_To_bw=","TSS=","houseKeeping=","c_proximal=","c_distal=","distal=","proximal="])




def usage():
	print( ' python3 idnetify_SOC.py -t input_file_for_call_peak -b bam_file -n Experiment name  -o OUTPUT_DIRECTORY  -g Effective_genome_size --g_file the_chrom_size_file\n\n -n Experiment name, which will be used to generate output file names  \n -p Number of processes to use \n --outdir OUTDIR \n --Broad_peak If --macs2_call_peak choose "off",please provide a broad peak bed file \n --footprint_file if choose off for --wellington_footprints,please provide the footprint file \n\n ###############Macs2 call peak \n\n -t input_file,the algin results \n -f Format of tag file,{AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}\n -g Effective genome size \n --pvalue the pvalue to filter peak,Default:5 \n --macs2_call_peak Whether using macs2 to call peak,can choose "on" or "off".Default:on \n\n ###############Footprint analysis \n\n -b ATAC_bam_file\n -A  Wether the data is ATAC data. \n --wellington_footprints Whether to execute footprint analysis,can choose "on" or "off".Default:on\n --Dnase_wig_To_bw If choose "on" for --wellington_footprints,this parameter can  choose "on" or "off" to decide whether to creat Dnase bw file.Default:on\n\n ###############Identify super ATAC broad peak \n\n --c_proximal the cut off of footrpitn plvalue and footprint number to identify super proximal ATAC peak \n --c_distal the cut off of footrpitn plvalue and footprint number to identify super distal ATAC peak \n --TSS the promoter file of gene \n --houseKeeping the promoter file of houseKeeping gene \n --proximal whether to identify proximal super peak,can choose "on" or "off".Default:on \n --distal whether to identify distal super peak,can choose "on" or "off".Default:on \n -h get help info')

delete_houseKeeping = "off"
macs2_call_peak_switch="on"
wellington_footprints_switch="on"
Dnase_wig_To_bw_switch="on"
proximal_identify= "on"
distal_identify= "on"
file_format = "AUTO"
ATAC_module =' ' 
Thread="1"
cut_off_proximal= "10"
cut_off_distal= "10"
peak_pvalue="5"
import os
default_path=os.getcwd()
output_path = default_path

for a,o in opts:
	if a in ('-t') :
		input_file = o
	if a in ('-f') :
		file_format = o
	if a in ('-g'):
		genome_size = o
	if a in ('-n'):
		output_name = o
	if a in ('-p'):
		Thread = o
	if a in ('-b'):
		input_bam_file = o
	if a in ('--outdir'):
		output_path = o
	if a in ('--g_file'):
		chrom_size_file = o
	if a in ('--Pvalue') :
		peak_pvalue = o
	if a in ('-A'):
		ATAC_module ="-A" 
	if a in ('-h'):
		usage()
		sys.exit()
	if a in ('--Broad_peak') :
		broad_peak_file_filter_sort = o
	if a in ('--Footprint_file') :
		footprint_file = o
	if a in ('--macs2_call_peak'):
		macs2_call_peak_switch= o
	if a in ('--wellington_footprints'):
		wellington_footprints_switch= o
	if a in ('--Dnase_wig_To_bw'):
		Dnase_wig_To_bw_switch = o
	if a in ('--TSS'):
		promoter_file = o
	if a in ('--houseKeeping'):
		houseKeeping_promoter_file = o 
		delete_houseKeeping = "on"
	if a in ('--c_proximal'):
		cut_off_proximal= o
	if a in ('--c_distal'):
		cut_off_distal= o
	if a in ('--proximal'):
		proximal_identify= o
	if a in ('--distal'):
		distal_identify= o

import subprocess
import datetime
import os.path
log_file1=open(output_path+"/Footprint_pipeline.log","a")
#############call broad peak
broad_peak_file=output_path+'/'+output_name+'_peaks.broadPeak'
broad_peak_file_filter=broad_peak_file+'_pfilter.bed'
if macs2_call_peak_switch== "on" :
    log_file2=open(output_path+"/macs2_call_peak.log","w")
    now_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print("###########################\n"+now_time+" Begain to call peak\n",file=log_file1)
    callpeak_command='macs2 callpeak  -t '+input_file +' -f '+file_format +' -g '+ genome_size+' -n '+ output_name+ ' --outdir '+ output_path+' -p 0.01 --nomodel --shift -75 --extsize 150 -B --broad  --keep-dup all '
    print(callpeak_command+"\n",file=log_file1)
    a1=subprocess.call(callpeak_command,shell=True,stderr=log_file2)
    now_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    if a1 ==0 :
        print(now_time+" Call peak complete\n",file=log_file1)
    else :
        print(now_time+" Call peak error,please read the macs2_call_peak.log \n",file=log_file1)
    subprocess.call("awk ' $1!~/\#/ && $8> "+ peak_pvalue +" {print}' " + broad_peak_file + ">" + broad_peak_file_filter,shell=True)
    broad_peak_file_filter_sort=broad_peak_file_filter+'.sort.bed'
    subprocess.call("bedtools sort -i " + broad_peak_file_filter + "|cut -f1,2,3 |awk '$3-$2 > 100 '>" + broad_peak_file_filter_sort,shell=True)
elif macs2_call_peak_switch== "off" :
    now_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print("###########################\n"+now_time+" Skip call peak\n",file=log_file1)
############footprint analysis
if wellington_footprints_switch=="on" :
    log_file3=open(output_path+"/wellington_footprints.log","w")
    footprint_out_path=output_path+"/wellington_footprints"
    if  os.path.exists(footprint_out_path) ==False :
        subprocess.call("mkdir "+ footprint_out_path ,shell=True)
    footprint_path=footprint_out_path+"/p_value_cutoffs/"+output_name
    now_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(now_time+" Begain to footprint analysis\n",file=log_file1)
    footprint_command="wellington_footprints.py -p " + Thread+" " + ATAC_module +" -o " + output_name +' ' +broad_peak_file_filter_sort +" " + input_bam_file +" " + footprint_out_path 
    print(footprint_command+"\n",file=log_file1)
    a2=subprocess.call(footprint_command,shell=True,stderr=log_file3)
    now_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    if a2 ==0 :
        print(now_time+" Footprint analysis complete\n",file=log_file1)
    else :
        print(now_time+" Footprint analysis error,please read the wellington_footprints.log \n",file=log_file1)
    subprocess.call("mv " + footprint_out_path +"/p\ value\ cutoffs "+footprint_out_path+ "/p_value_cutoffs ",shell=True)
    WellingtonFootprints_wig_file=footprint_out_path+"/"+output_name+".WellingtonFootprints.wig"
    subprocess.call("perl -e 'while(<>){$_=~s/start\=0/start\=1/ig if($_=~/[A-Z]/);print $_}' "+WellingtonFootprints_wig_file+" > "+footprint_out_path+"/"+output_name+".WellingtonFootprints_tmp.wig",shell=True)
    subprocess.call("wigToBigWig -clip "+footprint_out_path+"/"+output_name+".WellingtonFootprints_tmp.wig " + chrom_size_file +" "+footprint_out_path+"/"+output_name+".WellingtonFootprints.bw",shell=True)
    WellingtonFootprints_bw_file=footprint_out_path+"/"+output_name+".WellingtonFootprints.bw"
    WellingtonFootprints_bdg_file=footprint_out_path+"/"+output_name+".WellingtonFootprints.bdg"
    footprint_file=footprint_path+".WellingtonFootprints.-10.bed"
    subprocess.call("bigWigToBedGraph  "+ WellingtonFootprints_bw_file + " "+WellingtonFootprints_bdg_file,shell=True)
    subprocess.call("bedtools intersect -a  "+ footprint_file + " -b "+WellingtonFootprints_bdg_file + "-wa -wb > "+ footprint_path+"WellingtonFootprints_overlap_score.bed" ,shell=True)
    subprocess.call("bedtools sort -i  "+ footprint_path+"WellingtonFootprints_overlap_score.bed  >" + footprint_path+"WellingtonFootprints_overlap_score_sort.bed" ,shell=True)
    subprocess.call("bedtools merge -i  "+ footprint_path+"WellingtonFootprints_overlap_score_sort.bed -c 4,10,5  -o distinct,min,distinct >" + footprint_path+"WellingtonFootprints_score.bed" ,shell=True)
    footprint_score_file=footprint_path+"WellingtonFootprints_score.bed"
################cut_wig_TO_bw
    if Dnase_wig_To_bw_switch== "on" :
        bw_file_path=footprint_out_path+"/cut_wig"
        log_file4=open(output_path+"/cut_wig_TO_bw.log","w")
        if  os.path.exists(bw_file_path) ==False :
            subprocess.call("mkdir "+ bw_file_path ,shell=True)
        now_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print(now_time+" Begain to creat cut bw file\n",file=log_file1)
        creat_dnase_wig_command="dnase_wig_tracks.py "+ broad_peak_file_filter_sort +" " +input_bam_file+" "+ bw_file_path+"/"+output_name+"_danse_fr.wig "+" "+ bw_file_path+"/"+output_name+"_danse_rev.wig "+ATAC_module
        subprocess.call(creat_dnase_wig_command,shell=True,stderr=log_file4)
        merge_fr_rev_wig_command=" encode_ATAC-seq_wig_merge.pl -fr " +bw_file_path+"/"+output_name+"_danse_fr.wig -rev " +  bw_file_path+"/"+output_name+"_danse_rev.wig -out " +bw_file_path+"/"+output_name+"_danse_merge_fr_rev.wig"
        subprocess.call(merge_fr_rev_wig_command,shell=True,stderr=log_file4)
        wigToBigWig_command_1="wigToBigWig -clip "+bw_file_path+"/"+output_name+"_danse_merge_fr_rev.wig " + chrom_size_file +" "+bw_file_path+"/"+output_name+"_danse_merge_fr_rev.bw"
        a3=subprocess.call(wigToBigWig_command_1,shell=True,stderr=log_file4)
        print(creat_dnase_wig_command+"\n\n"+merge_fr_rev_wig_command+"\n\n"+wigToBigWig_command_1+"\n",file=log_file1)
        now_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        if a3 ==0 :
             print(now_time+" Creating of cut bw file complete\n",file=log_file1)
        else :
             print(now_time+" Creating of cut bw file error,please read the cut_wig_TO_bw.log \n",file=log_file1)
    elif Dnase_wig_To_bw_switch=="off" :
        now_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print(now_time+" Skip cut_wig_TO_bw analysis\n",file=log_file1)
elif  wellington_footprints_switch== "off" :
    now_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(now_time+" Skip footprint analysis\n",file=log_file1)

###################idnetify super broad ATAC peak
now_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print(now_time+" Begain to identify super broad peak\n",file=log_file1)
def identify_super_peak (ATAC_peak,footprint,cut_off,out_path) :
    subprocess.call(' Footprint_pvalue.R '+ footprint+" "+out_path+"/ATAC_peak_footprint_pvalue_stat.tiff "+out_path+"/footprint_pvalue_cut_off",shell=True)
    footprint_pvalue_cut_off=open(out_path+'/footprint_pvalue_cut_off',"r")
    data=footprint_pvalue_cut_off.read()
    data=data.strip("\n")
    data=data.split("\n")
    cut_off_dict1={}
    for line in  data :
        i=line.split("\t")
        cut_off_dict1[i[0]]=i[1]
    cut_off="top"+cut_off+"%"
    pvalue=cut_off_dict1[cut_off]
    footprint_pvalue_filter=out_path+"/"+output_name+".WellingtonFootprints."+pvalue+".bed"
    subprocess.check_call("awk '$5< " +pvalue + "{print $0}' " + footprint +" > "+ footprint_pvalue_filter ,shell=True)
    peak_overlap_footprint=open(out_path+"/"+output_name+"_broadPeak_overlap_footprints.bed","w")
    subprocess.check_call(["bedtools","intersect","-a",ATAC_peak,"-b",footprint_pvalue_filter,"-wa","-wb"],stdout=peak_overlap_footprint,shell=False)
    subprocess.call("awk '{print $1"+ '"-" $2'+ '"-"'+ "$3}' " +out_path+"/"+ output_name+"_broadPeak_overlap_footprints.bed"+ " | sort |uniq -c  |awk '{print $2 "+'"\t"'+" $1}' > " +out_path+"/"+output_name+"_broadPeak_footprint_number",shell=True)
    subprocess.call(' Footprint_number.R '+ out_path+"/"+output_name+"_broadPeak_footprint_number "+out_path+"/ATAC_peak_footprint_number_stat.tiff "+ out_path+"/footprint_number_cut_off",shell=True)
    footprint_number_cut_off=open(out_path+'/footprint_number_cut_off',"r")
    data1=footprint_number_cut_off.read()
    data1=data1.strip("\n")
    data1=data1.split("\n")
    cut_off_dict2={}
    for line in  data1 :
        i=line.split("\t")
        cut_off_dict2[i[0]]=i[1]
    number=cut_off_dict2[cut_off]
    subprocess.check_call("awk '$2>= "+number+' {print $1 '+ '"\t"'+ "$2}' "+ out_path+"/"+output_name+"_broadPeak_footprint_number | sed 's/-/"+r'\t'+"/g' > "+ out_path+ "/"+ output_name +"_super_ATAC_broad_peak.bed" ,shell=True)
    subprocess.check_call("awk '$2< "+number+' {print $1 '+ '"\t"'+ "$2}' "+ out_path+"/"+output_name+"_broadPeak_footprint_number | sed 's/-/"+r'\t'+"/g' > "+ out_path+ "/"+ output_name +"_typical_ATAC_broad_peak.bed" ,shell=True)
    subprocess.check_call("bedtools intersect -a "+ ATAC_peak+" -b " + footprint_pvalue_filter+ " -wa -v |awk '{print $1 "+ ' "\t" $2'+'"\t" $3 "\t"'+" 0} ' |cat >> " + out_path+ "/"+ output_name +"_typical_ATAC_broad_peak.bed",shell=True)
    return
####distal
if distal_identify == "on" : 
    distal_out_path=output_path+"/distal"
    distal_ATAC_peak=output_path+"/" +output_name +"_distal_ATAC_broad_peak.bed"
    distal_footprint=output_path+"/" + output_name +"_distal.WellingtonFootprints.10.bed"
    if  os.path.exists(distal_out_path) ==False :
        subprocess.call("mkdir "+ distal_out_path ,shell=True)
    subprocess.check_call("bedtools intersect -a "+ broad_peak_file_filter_sort+" -b " + promoter_file + " -wa -v > " + distal_ATAC_peak,shell=True)
    subprocess.check_call("bedtools intersect -a "+footprint_score_file +" -b "  + distal_ATAC_peak + " -wa -u > " +distal_footprint ,shell=True)
    identify_super_peak(distal_ATAC_peak,distal_footprint,cut_off_distal,distal_out_path)
####proximal
if proximal_identify == "on" :
    proximal_out_path=output_path+"/proximal"
    proximal_ATAC_peak=output_path+"/" +output_name +"_proximal_ATAC_broad_peak.bed"
    proximal_ATAC_peak_non_HK=output_path+"/" +output_name +"_proximal_ATAC_broad_peak_non_HK_gene.bed"
    proximal_footprint=output_path+"/" + output_name +"_proximal.WellingtonFootprints.10.bed"
    proximal_footprint_non_HK=output_path+"/" + output_name +"_proximal_non_HK_gene.WellingtonFootprints.-10.bed"
    if  os.path.exists(proximal_out_path) ==False :
        subprocess.call("mkdir "+ proximal_out_path ,shell=True)
    if delete_houseKeeping == "off" :
        subprocess.check_call("bedtools intersect -a "+ broad_peak_file_filter_sort+" -b " + promoter_file + " -wa -u > " + proximal_ATAC_peak,shell=True)
        subprocess.check_call("bedtools intersect -a "+footprint_score_file +" -b "  + proximal_ATAC_peak + " -wa -u > " +proximal_footprint ,shell=True)
        identify_super_peak(proximal_ATAC_peak,proximal_footprint,cut_off_proximal,proximal_out_path)
    elif delete_houseKeeping == "on" :
        subprocess.check_call("bedtools intersect -a "+ broad_peak_file_filter_sort+" -b " + promoter_file + " -wa -u > " + proximal_ATAC_peak,shell=True)
        subprocess.check_call("bedtools intersect -a "+ proximal_ATAC_peak+" -b " + houseKeeping_promoter_file + " -wa -v > " + proximal_ATAC_peak_non_HK,shell=True)
        subprocess.check_call("bedtools intersect -a "+footprint_score_file +" -b "  + proximal_ATAC_peak + " -wa -u > " +proximal_footprint ,shell=True)
        subprocess.check_call("bedtools intersect -a "+footprint_score_file +" -b "  + proximal_ATAC_peak_non_HK + " -wa -u > " +proximal_footprint_non_HK,shell=True)
        identify_super_peak(proximal_ATAC_peak_non_HK,proximal_footprint_non_HK,cut_off_proximal,proximal_out_path)
now_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print(now_time+" Identifying of super broad peak complete\n",file=log_file1)
