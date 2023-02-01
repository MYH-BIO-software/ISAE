# ISAE: Identify super accessibility elements  

## 1.Introduction  
The software was used to identify super accessibility elements (SAEs) using ATAC-seq data, which were open chromatin regions in promoters or enhancers that are strongly bound simultaneously by multiple lineage-specific TFs. In addition to this, the software can also associate SAEs with gene expression and identify differential SAEs between different tissues/cell stages and types considering the change both in ATAC signal and TF footprint binding.   

## 2.Dependencies  
python3  
R > v3.4  
perl > v5  
bedtools > v2  
samtools > v1  
macs2 > v2  
pyDNase v0.030  
HOMER > v4  
wigToBigWig  
bigWigToBedGraph  

## 3.Usage  
### 1>Idnetify_SAE.py  
python3 Idnetify_SAE.py [options] -t Input_file_for_call_peak -b Bam_file -n Experiment_name  -o Output_directory  -g Effective_genome_size --g_file the_chrom_size_file  

 -n Experiment name, which will be used to generate output file names  
 -p Number of processes to use  
 --outdir OUTDIR  
 --Broad_peak If --macs2_call_peak choose "off",please provide a broad peak bed file  
 --footprint_file if choose off for --wellington_footprints,please provide the footprint file  

 ###############Macs2 call peak  

 -t input_file,the algin results  
 -f Format of tag file,{AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}  
 -g Effective genome size  
 --pvalue the pvalue to filter peak,Default:5  
 --macs2_call_peak Whether using macs2 to call peak,can choose "on" or "off".Default:on  

 ###############Footprint analysis  

 -b ATAC_bam_file  
 -A  Wether the data is ATAC data.  
 --wellington_footprints Whether to execute footprint analysis,can choose "on" or "off".Default:on  
 --Dnase_wig_To_bw If choose "on" for --wellington_footprints,this parameter can  choose "on" or "off" to decide whether to creat Dnase bw file.Default:on  

 ###############Identify super ATAC broad peak  

 --c_proximal the cut off of footrpitn plvalue and footprint number to identify super proximal ATAC peak  
 --c_distal the cut off of footrpitn plvalue and footprint number to identify super distal ATAC peak  
 --TSS the promoter file of gene  
 --houseKeeping the promoter file of houseKeeping gene  
 --proximal whether to identify proximal super peak,can choose "on" or "off".Default:on  
 --distal whether to identify distal super peak,can choose "on" or "off".Default:on  
 -h Print this help menu  

### 2>link_ATAC_peak_with_gene.py  
python3 link_ATAC_peak_with_gene.py -s super peak bed file -t typical peak bed file -g  gene tss bed file  -G gff annotation file  -k houseKeepgene_ID_list -h Print this help menu  

### 3>DifferentialSAE.py  
python3 DifferentialSAE.py  -pt  Treatment ATAC broad peak bed file -pc bam_file  Control ATAC broad bed file -bt  Treatment bam file  -bc  Control bam file  -st  Treatment SAE  bed file  -sc  Control SAE  bed file -nt Treatment sample name -nc Control sample name -o output path  

-pt ATAC broad peak bed file of treatment group  
-pc ATAC broad peak bed file of control group  
-bt ATAC bam file of treatment group  
-bc ATAC bam file of control group  
-st SAE  bed file of treatment group -st SAE  bed file of control group  
-nt sample name of treatment group  
-nc sample name of control group  
-o Path to output folder  
-h Print this help menu  
