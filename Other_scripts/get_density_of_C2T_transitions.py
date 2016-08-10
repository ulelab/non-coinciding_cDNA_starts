'''
Created on June 25, 2016
 
@author: Nejc Haberman
 
description:
Script will recive a special SAM file format (read bellow) and return a density array of positions where C to T transitions are located across all cDNAs
For more details about SAM format: https://samtools.github.io/hts-specs/SAMv1.pdf
 
SAM format preparation: 
- samtools are needed http://samtools.sourceforge.net/
- we used samtools valmd -e option, which replace each nucleotide of the sequence that matches with genomic sequence by replacing it with '=' and all mutations are left
    - example: samtools calmd -e -u example.bam genomic_sequence.fasta > example_with_transitions.sam
'''

import sys

def get_deletions_coverage (fin_sam):
    fin = open(fin_sam, "rt")
    line = fin.readline()
    deletion_sum = 0    #number of all deletions
    read_count = 0
    coverage = 150 * [0] #max cDNA length is 150 nt
    deletion_reads_count = 0
    while line[0] == '@':    #ignore header of .SAM
        line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
        strand = col[1] # 0 is plus strand, 16 is minus
        seq = col[9]    #cDNA sequnce that matches the genomic sequence is now replaced with '=' and only mutations are left
        read_count += 1
        if strand == "0":   #plus strand
            trans_pos = [n for n in xrange(len(seq)) if seq.find('C', n) == n]  #get all positions of C transitions
            for i in range(0,trans_pos.__len__()): #if we have a multiple C transitions they can be one after another
                pos = trans_pos[i]
                coverage[pos] += 1
        if strand == "16":  #minus strand
            seq = seq[::-1] #because of the minust strand we turned around the sequence
            trans_pos = [n for n in xrange(len(seq)) if seq.find('G', n) == n]  #get all positions of transitions to G
            for i in range(0,trans_pos.__len__()): #if we have a multiple G transitions they can be one after another
                pos = trans_pos[i]
                coverage[pos] += 1
        if trans_pos != -1:
            deletion_reads_count += 1
    line = fin.readline()
        
    print "\nDensity of C transitions across all cDNAs:\n"
    print str(coverage).replace('[','').replace(']','').replace("'","")
    print str(deletion_sum)
    print "\nNumber of cDNAs with C transitions:\t" + str(deletion_reads_count)
    print "\nNumber of all cDNAs:\t" + str(read_count)
    fin.close()

# main
if sys.argv.__len__() == 2:
    fin_sam = sys.argv[1]
    get_deletions_coverage(fin_sam)
else:
    print "you need 1 argument to run the script in preprocessed SAM format"
    quit()

