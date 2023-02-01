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

Most of these dependencies can be installed at once using conda: http://conda.pydata.org/miniconda.html  
## 3.Installation  
```
git clone https://github.com/MingyangHu/ISAE.git
export PATH=$PATH:/path/to/SAE/bin
```

## 4.Usage  
### 1>Idnetify_SAE.py  
The script was used to call ATAC broad peak, identify footprints in ATAC peak and super accessibility elements.    
```
Idnetify_SAE.py [options] -t Input_file_for_call_peak -b Bam_file -n Experiment_name  -o Output_directory  -g Effective_genome_size --g_file the_chrom_size_file  
```
* Explanation of parameters  
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
--Dnase_wig_To_bw If choose "on" for  
--wellington_footprints,this parameter can  choose "on" or "off" to decide whether to creat Dnase bw file.Default:on  
###############Identify SAE  
--c_proximal the cut off of footrpitn plvalue and footprint number to identify super proximal ATAC peak  
--c_distal the cut off of footrpitn plvalue and footprint number to identify super distal ATAC peak  
--TSS the promoter file of gene  
--houseKeeping the promoter file of houseKeeping gene  
--proximal whether to identify proximal super peak,can choose "on" or "off".Default:on  
--distal whether to identify distal super peak,can choose "on" or "off".Default:on  
-h Print this help menu  

### 2>link_ATAC_peak_with_gene.py  
The script was used to assign super and typical accessibility elements to each expressed gene and compare the gene expression between the two group genes.    
```
link_ATAC_peak_with_gene.py -s super peak bed file -t typical peak bed file -g  gene tss bed file  -G gff annotation file  -k houseKeepgene_ID_list -h Print this help menu  
```
### 3>DifferentialSAE.py  
The script was used to identify differential SAEs between different tissues/cell stages and types considering the change both in ATAC signal and TF footprint binding.     
```
DifferentialSAE.py  -pt  Treatment_ATAC_broad_peak_bed_file -pc   Control_ATAC_broad_bed_file -bt  Treatment_bam_file  -bc  Control_bam_file  -st  Treatment_SAE_bed_file  -sc  Control_SAE_bed_file -nt Treatment_sample_name -nc Control_sample_name -o output_path  
```
* Explanation of parameters  
-pt ATAC broad peak bed file of treatment group  
-pc ATAC broad peak bed file of control group  
-bt ATAC bam file of treatment group  
-bc ATAC bam file of control group  
-st SAE  bed file of treatment group -st SAE  bed file of control group  
-nt sample name of treatment group  
-nc sample name of control group  
-o Path to output folder  
-h Print this help menu  
