#!/bin/bash -l
#$ -l h_vmem=16G
#$ -l tmem=16G
#$ -l h_rt=32:0:0
#$ -j y
#$ -S /bin/bash

# list of all tools
'''
PYTHONPATH=/home/programs/Python-2.7.5/bin
PERL5LIB=/home/programs/ActivePerl-5.16.3.1603-x86_64-linux-glibc-2.3.5-296746/perl/bin
export PATH=$PYTHONPATH:$PATH
export PATH=$PERL5LIB:$PATH
export PATH=/home/programs/tophat-2.0.9.Linux_x86_64:$PATH
export PATH=/home/programs/bowtie2-2.1.0:$PATH
export PATH=/home/programs/samtools-0.1.19:$PATH
export PATH=/home/programs/fastx_toolkit0.0.13:$PATH
export PATH=/home/programs/bedtools-2.17.0/bin:$PATH
'''
data=$1
path=`pwd -P`

gunzip ${path}${data}.gz

# clip the adapter and discard clipped sequences and discard the sequences that are shorter then 17 nt; complete cDNAs contains the adapter sequence and incomplete does not
fastx_clipper -Q 33 -a AGATCGGAAG -C -n -l 17 -i  ${path}${data} -o  ${path}${data}-incomplete.fq
fastx_clipper -Q 33 -a AGATCGGAAG -c -n -l 17 -i  ${path}${data} -o  ${path}${data}-complete.fq

# fastq to fasta
fastq_to_fasta -Q 33 -n -i ${path}${data}-incomplete.fq -o ${path}${data}-incomplete.fa
fastq_to_fasta -Q 33 -n -i ${path}${data}-complete.fq -o ${path}${data}-complete.fa
rm ${path}${data}-incomplete.fq
rm ${path}${data}-complete.fq

# map to hg19
bowtie2-align -x ~/bowtie-indexes/hg19/hg19 -f ${path}${data}-incomplete.fa -S ${path}${data}-incomplete.sam
bowtie2-align -x ~/bowtie-indexes/hg19/hg19 -f ${path}${data}-complete.fa -S ${path}${data}-complete.sam
rm ${path}${data}-incomplete-barcodes.fa
rm ${path}${data}-complete-barcodes.fa

# filter reads with more then 2 mismatches
samtools view -Sh ${path}${data}-incomplete.sam | grep -e "^@" -e "XM:i:[012][^0-9]" > ${path}${data}-incomplete-2mis.sam
samtools view -Sh ${path}${data}-complete.sam | grep -e "^@" -e "XM:i:[012][^0-9]" > ${path}${data}-complete-2mis.sam

# SAM to BAM
samtools view -hSb ${path}${data}-incomplete-2mis.sam > ${path}${data}-incomplete-2mis.bam
samtools view -hSb ${path}${data}-complete-2mis.sam > ${path}${data}-complete-2mis.bam
rm ${path}${data}-incomplete-2mis.sam
rm ${path}${data}-complete-2mis.sam

# convert bam to bed
bedtools bamtobed -i ${path}${data}-incomplete-2mis.bam > ${path}${data}-long-2mis.bed
bedtools bamtobed -i ${path}${data}-complete-2mis.bam > ${path}${data}-short-2mis.bed

# compress the original data
gzip ${path}${data}

