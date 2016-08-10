#!/bin/bash -l
# This is the main script for finding cDNA start and cDNA end peaks for eIFA3 data
# Input:	mapped complete cDNAs (cDNAs that contains the adapter sequence) to transcriptome in BED format
# Output:	a tab delimited text file with three columns: junction(ID:position), number of cDNAs, a peak position ralative to junction

cDNAs=$1

# filter exons smaller then 100 nts
python getRNAmapPositions-both-directions-filter100_100.py ${cDNAs} ${cDNAs}-filtered100.bed

# select cDNAs from top 1000 transcripts
Rscript select_top1000_transcripts.R ${cDNAs} ${cDNAs}-selected.bed

# get peaks
Rscript get_peaks-start_peaks.R ${cDNAs}-selected.bed ${cDNAs}-selected-start_peaks.tab
Rscript get_peaks-end_peaks.R ${cDNAs}-selected.bed ${cDNAs}-selected-end_peaks.tab

# plot a density figure of cDNA starts and cDNA ends around peaks
Rscript draw_a_density_map_around_peak_reference.R ${cDNAs}-selected-start_peaks.tab ${cDNAs}
Rscript draw_a_density_map_around_peak_reference.R ${cDNAs}-selected-end_peaks.tab ${cDNAs}
