#!/bin/bash -l
#$ -l h_vmem=16G
#$ -l tmem=16G
#$ -l h_rt=32:0:0
#$ -j y
#$ -S /bin/bash

PYTHONPATH=/home/programs/Python-2.7.5/bin
PERL5LIB=/home/programs/ActivePerl-5.16.3.1603-x86_64-linux-glibc-2.3.5-296746/perl/bin
export PATH=$PYTHONPATH:$PATH
export PATH=$PERL5LIB:$PATH
export PATH=/home/programs/tophat-2.0.9.Linux_x86_64:$PATH
export PATH=/home/programs/bowtie2-2.1.0:$PATH
export PATH=/home/programs/samtools-0.1.19:$PATH
export PATH=/home/programs/fastx_toolkit0.0.13:$PATH
export PATH=/home/programs/bedtools-2.17.0/bin:$PATH

data=$1
path=`pwd -P`

gunzip ${path}${data}.gz

# remove adapter seq and keep all reads that are longer then 26 nts which 17 nts (17nts seq + 5nts random barcode + 4nts experimental barcode); complete cDNAs contains the adapter sequence and incomplete does not
fastx_clipper -Q 33 -C -n -l 26 -a AGATCGGAAG -i ${path}${data} -o ${path}${data}-incomplete.fq
fastx_clipper -Q 33 -c -n -l 26 -a AGATCGGAAG -i ${path}${data} -o ${path}${data}-complete.fq

# convert fastq to fasta
fastq_to_fasta -Q 33 -n -i ${path}${data}-incomplete.fq -o ${path}${data}-incomplete.fa
fastq_to_fasta -Q 33 -n -i ${path}${data}-complete.fq -o ${path}${data}-complete.fa
rm ${path}${data}-incomplete.fq
rm ${path}${data}-complete.fq

# swap random barcodes to headers of fasta file
python ${path}swap_barcode_to_header.py ${path}${data}-incomplete.fa ${path}${data}-incomplete-barcodes.fa
python ${path}swap_barcode_to_header.py ${path}${data}-complete.fa ${path}${data}-complete-barcodes.fa
rm ${path}${data}-incomplete.fa
rm ${path}${data}-complete.fa

# map to transcripts
# important: this is a custom made bowtie index build from Ensemble Biomart by selecting fasta from only coding sequences (CDS) and saving exon positions in the header of the fasta
bowtie2-align -x /home/bowtie-indexes/Human-GRCh38.p2-CDS-transcripts/Human-GRCh38.p2-CDS-transcripts -f ${path}${data}-incomplete-barcodes.fa -S ${path}${data}-incomplete.sam
bowtie2-align -x /home/bowtie-indexes/Human-GRCh38.p2-CDS-transcripts/Human-GRCh38.p2-CDS-transcripts -f ${path}${data}-complete-barcodes.fa -S ${path}${data}-complete.sam
rm ${path}${data}-incomplete-barcodes.fa
rm ${path}${data}-complete-barcodes.fa

# filter out all reads with more then 2 mismatches
samtools view -Sh ${path}${data}-incomplete.sam | grep -e "^@" -e "XM:i:[012][^0-9]" > ${path}${data}-incomplete-2mis.sam
samtools view -Sh ${path}${data}-complete.sam | grep -e "^@" -e "XM:i:[012][^0-9]" > ${path}${data}-complete-2mis.sam
rm ${path}${data}-incomplete.sam
rm ${path}${data}-complete.sam

# SAM to BAM
samtools view -hSb ${path}${data}-incomplete-2mis.sam > ${path}${data}-incomplete-2mis.bam
samtools view -hSb ${path}${data}-complete-2mis.sam > ${path}${data}-complete-2mis.bam
rm ${path}${data}-incomplete-2mis.sam
rm ${path}${data}-complete-2mis.sam

# convert bam to bed
bedtools bamtobed -i ${path}${data}-incomplete-2mis.bam > ${path}${data}-incomplete-2mis.bed
bedtools bamtobed -i ${path}${data}-complete-2mis.bam > ${path}${data}-complete-2mis.bed

# remove duplicates
cat ${path}${data}-incomplete-2mis.bed | sort -k1,1 -k2,2n -k5,5 | uniq > ${path}${data}-long-2mis-uniq.bed
cat ${path}${data}-complete-2mis.bed | sort -k1,1 -k2,2n -k5,5 | uniq > ${path}${data}-short-2mis-uniq.bed

# compress the original data
gzip ${path}${data}



