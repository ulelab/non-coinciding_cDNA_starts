#!/bin/bash -l
# Clusters must be in BED foramt and kmers separated by new lines. The fasta file was downlaoded from ucsc website for hg19 genome.
# Tested on bedtools version 2.22.1 

# input
clusters=$1
kmers=$2
name=$3
fasta=../ucsc.hg19.fasta

# sorted clusters by their length
Rscript get_cluster_length.R ${clusters} ${clusters}-length.bed

# get flanked surrounding positions 
python flank_clusters.py ${clusters}-length.bed ${clusters}-flanked50_200.bed 50 200

# get fasta 
bedtools getfasta -s -name -fi /media/skgthab/storage/UCL/genome-data/ucsc.hg19.fasta -bed ${clusters}-flanked50_200.bed -fo ${clusters}-flanked50_200.fasta

# get kmer enrichment in csv HeatMap format for lenght grouped clusters
python kmer_coverage-grouped_clusters.py ${clusters}-flanked50_200.fasta ${kmers} ${name}.csv

# draw a HeatMap
Rscript draw-HeatMap.R ${name}.csv ${name}.pdf ${name}


