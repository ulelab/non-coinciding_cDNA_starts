#!/bin/bash -l
# 
# To identify cross-linked clusters we will first have to define significantly cross linked binding sites by using iCount tools:
# https://github.com/tomazc/iCount
# 
# Tested on bedtools version 2.22.1 

# we are using a virtual environment for iCount: http://icount.readthedocs.io/en/latest/installation.html#installing-from-source
source ../iCount/bin/activate

genome=gene.gtf.gz	#gene positions across genome to normalise peaks on the single transcript level
data=$1	#mapped cDNAs in BED format
hw=3	#half window size

mkdir ${data}${hw}nt.tmp #temp folder
export ICOUNT_TMP_ROOT=${data}${hw}nt.tmp #export temp folder for iCount

# get cross-link positions from mapped cDNAs that are in BED format
python BEDtoXlink.py ${data} ${data}.xlink.bed
sort -k1,1 -k2,2n -k6,6 ${data}.xlink.bed > ${data}${hw}nt.tmp/xlink-sorted.tmp
python BEDsum-iCount.py ${data}${hw}nt.tmp/xlink-sorted.tmp > ${data}.xlink.sum.bed

# find significant cross-link peaks and merge them into clusters
iCount peaks --half_window ${hw} --fdr 0.05 --perms 100 ${genome} ${data}.xlink.sum.bed ${data}.xlink.sum.${hw}nt.hw-peaks.tab --scores ${data}.xlink.sum.${hw}nt.hw-scores.tsv
cat ${data}.xlink.sum.${hw}nt.hw-peaks.tab | awk '{print $1 "\t" $2 "\t" $3 "\t.\t.\t" $8}' > ${data}.xlink.sum.${hw}nt.hw-peaks.bed	#convert it to bed format
bedtools merge -i ${data}.xlink.sum.${hw}nt.hw-peaks.bed -s -d 3 -c 5,5,6 -o distinct,count,distinct > ${data}.xlink.sum.${hw}nt.hw-peaks.3nt.clusters.bed

rm -r ${data}${hw}nt.tmp
