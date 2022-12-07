# PCA-pipline
PCA pipeline for whole genome data
read me, how I did it, scripts, figures, metadata


# Table of contents



# Obtaining the data
Most of the raw reads were obtained using illumina sequencing data downloaded from the NCBI SRA database. Reads in this analysis that were not downloaded from a database were prepped by the Uehling lab at OSU and sequenced by the UO sequencing core using illumina. The reference genome used was obtained from the NCBI SRA database https://www.ncbi.nlm.nih.gov/Traces/wgs/SMRR01?display=contigs&page=1. 

# Aligning short reads to reference genome
BWA manual / web documentation
http://bio-bwa.sourceforge.net/bwa.shtml

bwa mem -t 10 $REFERENCEGENOME ${read1} ${read2} > output.sam 


samtools view -bS output.sam > output.bam
samtools sort -@ output.bam -o 2ndoutput.bam
