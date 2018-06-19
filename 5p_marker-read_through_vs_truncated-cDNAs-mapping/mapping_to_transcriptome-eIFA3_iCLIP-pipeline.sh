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

# clip the adapter and discard clipped sequences and discard the sequences that are shorter then 17 nt + 5 random barcode + 4 experimental barcode; complete cDNAs contains the adapter sequence and incomplete does not
fastx_clipper -Q 33 -a AGATCGGAAG -n -l 26 -i  ${path}${data} -o  ${path}${data}-trimmed.fq

# fastq to fasta with random barcodes to headers and separation of read-through and truncated reads 
python ${path}FASTQ2FASTA-swap_barcode_to_header-readThrough.py ${path}${data}-trimmed.fq ${path}${data}-trimmed

# map to hg19
bowtie2-align -x /home/bowtie-indexes/Human-GRCh38.p2-CDS-transcripts/Human-GRCh38.p2-CDS-transcripts -f ${path}${data}-trimmed-readThrough.fa -S ${path}${data}-trimmed-readThrough.sam
bowtie2-align -x /home/bowtie-indexes/Human-GRCh38.p2-CDS-transcripts/Human-GRCh38.p2-CDS-transcripts -f ${path}${data}-trimmed-truncated.fa -S ${path}${data}-trimmed-truncated.sam
rm ${path}${data}-trimmed-readThrough.fa
rm ${path}${data}-trimmed-truncated.fa

# filter reads with more then 2 mismatches
samtools view -Sh ${path}${data}-trimmed-readThrough.sam | grep -e "^@" -e "XM:i:[012][^0-9]" > ${path}${data}-trimmed-readThrough-2mis.sam
samtools view -Sh ${path}${data}-trimmed-truncated.sam | grep -e "^@" -e "XM:i:[012][^0-9]" > ${path}${data}-trimmed-truncated-2mis.sam

# SAM to BAM
samtools view -hSb ${path}${data}-trimmed-readThrough-2mis.sam > ${path}${data}-trimmed-readThrough-2mis.bam
samtools view -hSb ${path}${data}-trimmed-truncated-2mis.sam > ${path}${data}-trimmed-truncated-2mis.bam
rm ${path}${data}-trimmed-readThrough-2mis.sam
rm ${path}${data}-trimmed-truncated-2mis.sam

# convert bam to bed
bedtools bamtobed -i ${path}${data}-trimmed-readThrough-2mis.bam > ${path}${data}-trimmed-readThrough-2mis.bed
bedtools bamtobed -i ${path}${data}-trimmed-truncated-2mis.bam > ${path}${data}-trimmed-truncated-2mis.bed

# remove duplicates
cat ${path}${data}-trimmed-readThrough-2mis.bed | sort -k1,1 -k2,2n -k5,5 | uniq > ${path}${data}-trimmed-readThrough-2mis-uniq.bed
cat ${path}${data}-trimmed-truncated-2mis.bed | sort -k1,1 -k2,2n -k5,5 | uniq > ${path}${data}-trimmed-truncated-2mis-uniq.bed
rm ${path}${data}-trimmed-readThrough-2mis-uniq.bed
rm ${path}${data}-trimmed-truncated-2mis-uniq.bed

# compress the original data
gzip ${path}${data}



