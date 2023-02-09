#!/usr/bin/python
import getopt, sys

opts, args = getopt.getopt(sys.argv[1:], "t:c:i:o:h",["nt=","nc="])

def usage():
	print(sys.argv[0] + ' -t treatment_footprint_number  -c control_footprint_number -i merge_peak  -o output_file  --nt Treatment sample name --nc  Control sample name ')
	print(sys.argv[0] + ' -h #get help info')


for a,o in opts:
	if a in ('-i'):
		input_file = o
	elif a in ('-t'):
		treatment_footprint_number = o
	elif a in ('-c'):
		control_footprint_number = o
	elif a in ('-o'):
		output_file = o
	elif a in ('--nc'):
		control_name = o
	elif a in ('--nt'):
		treatment_name = o
	elif a in ('-h'):
		usage()
		sys.exit()

file1=open(treatment_footprint_number,"r")
data1=file1.readlines()
treatment_footprint_Num={}
for line in data1 :
    line=line.strip("\n")
    i=line.split("\t")
    treatment_footprint_Num[i[0]]=i[1]

file2=open(control_footprint_number,"r")
data2=file2.readlines()
control_footprint_Num={}
for line in data2 :
    line=line.strip("\n")
    i=line.split("\t")
    control_footprint_Num[i[0]]=i[1]

file3=open(input_file,"r")
file4=open(output_file,"w")
data3=file3.readlines()
for line in data3 :
    line=line.strip("\n")
    i=line.split("\t")
    peak_list=i[3].split(",")
    peak_num=len(peak_list)
    if int(peak_num) == 1 :
        if control_name in i[3] :
            print(i[3]+"\t"+control_footprint_Num[i[3]]+"\t"+"0"+"\t"+ str(0-int(control_footprint_Num[i[3]])) ,file=file4)
        elif treatment_name in i[3] :
            print(i[3]+"\t"+"0"+"\t"+treatment_footprint_Num[i[3]]+"\t"+ str(int(treatment_footprint_Num[i[3]])-0),file=file4)
    elif int(peak_num) == 2 :
        a=i[3].split(",")
        if control_name in  a[1] : 
            print(i[3]+"\t"+control_footprint_Num[a[1]]+"\t"+treatment_footprint_Num[a[0]]+"\t"+str(int(treatment_footprint_Num[a[0]])-int(control_footprint_Num[a[1]])),file=file4)
        else :
            print(i[3]+"\t"+control_footprint_Num[a[0]]+"\t"+treatment_footprint_Num[a[1]]+"\t"+str(int(treatment_footprint_Num[a[1]])-int(control_footprint_Num[a[0]])),file=file4)
