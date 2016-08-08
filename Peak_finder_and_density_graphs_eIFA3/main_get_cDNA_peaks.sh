#!/bin/bash -l
# This is the main script for finding cDNA start and cDNA end peaks for eIFA3 data
# Input:	mapped cDNAs to transcriptome in BED format
# Output:	a tab delimited text file with three columns: junction(ID:position), number of cDNAs, a peak position ralative to junction

reads1=$1

# filter exons smaller then 100 nts
python getRNAmapPositions-both-directions-filter100_100.py ${reads1} ${reads1}-filtered100.bed

# select cDNAs from top 1000 transcripts
Rscript select_topN_transcripts.R ${reads1} ${reads1}-selected.bed

# get peaks
Rscript get_peaks-start_peaks.R ${reads1}-selected.bed ${reads1}-selected-start_peaks.tab
Rscript get_peaks-end_peaks.R ${reads1}-selected.bed ${reads1}-selected-end_peaks.tab

# plot a density figure of cDNA starts and cDNA ends around peaks
Rscript draw_a_density_map_around_peak_reference.R ${reads1}-selected-start_peaks.tab
Rscript draw_a_density_map_around_peak_reference.R ${reads1}-selected-end_peaks.tab
