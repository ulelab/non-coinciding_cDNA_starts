#!/bin/bash -l

#This is the main script for Normaisation and density graph of PTBP1. 
# Input: clusters of your interest, mapped reads in size of your interest

reference=$1
cDNAs=$2
name="Normalisation_and_density_graph-PTBP1"

#get cluster start and end position
pytgon getStartAndEnd-BED.py ${reference} ${reference}

#flank reference by 120 nt
python flankBEDpositionsCustom.py ${reference}-Start.bed ${reference}-Start-flanked120.bed 120 120
python flankBEDpositionsCustom.py ${reference}-End.bed ${reference}-End-flanked120.bed 120 120

#get start and end positions from dDNAs (reads)
python getStartAndEnd-BED.py ${cDNAs} ${cDNAs}

#get read-starts around cluster start and end position
bedtools intersect -s -a ${cDNAs}-Start.bed -b ${reference}-Start-flanked120.bed -wb > 01.tmp
bedtools intersect -s -a ${cDNAs}-Start.bed -b ${reference}-Start-flanked120.bed -wb > 02.tmp
bedtools intersect -s -a ${cDNAs}-End.bed -b ${reference}-End-flanked120.bed -wb > 03.tmp
bedtools intersect -s -a ${cDNAs}-End.bed -b ${reference}-End-flanked120.bed -wb > 04.tmp

#get original BP position back
cat 01.tmp | tr ':' '\t' > 01.tab.tmp
cat 02.tmp | tr ':' '\t' > 02.tab.tmp
cat 03.tmp | tr ':' '\t' > 03.tab.tmp
cat 04.tmp | tr ':' '\t' > 04.tab.tmp

#plot maps
Rscript Normalisation_and_density_graph-PTBP1.R 01.tab.tmp 02.tab.tmp 03.tab.tmp 04.tab.tmp ${name}

#clean
rm *.tmp
rm ${reference}-Start-flanked120.bed
rm ${reference}-End-flanked120.bed
rm ${cDNAs}-Start.bed
rm ${cDNAs}-End.bed
