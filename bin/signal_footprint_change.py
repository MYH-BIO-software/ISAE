#!/usr/bin/python
import getopt, sys

opts, args = getopt.getopt(sys.argv[1:], "s:f:o:h")

def usage():
	print(sys.argv[0] + ' -s peak singal change file  -f peak footprint change file   -o output_file  ')
	print(sys.argv[0] + ' -h #get help info')


for a,o in opts:
	if a in ('-s'):
		signal_chnage = o
	elif a in ('-f'):
		footprint_chnage = o
	elif a in ('-o'):
		output_file = o
	elif a in ('-h'):
		usage()
		sys.exit()


file1=open(signal_chnage,"r")
data1=file1.readlines()
signal_change_dict={}
for line in data1 :
    line=line.strip("\n")
    i=line.split("\t")
    signal_change_dict[i[0]]=i[1]+"\t"+i[2]+"\t"+i[3]+"\t"+i[0]+"\t"+i[4]+"\t"+i[5]

file2=open(footprint_chnage,"r")
data2=file2.readlines()
footprint_change_dict={}
for line in data2 :
    line=line.strip("\n")
    i=line.split("\t")
    footprint_change_dict[i[0]]=i[3]

file3=open(output_file,"w")

for i in signal_change_dict :
    if i in footprint_change_dict :
        print(signal_change_dict[i]+"\t"+footprint_change_dict[i],file=file3)
