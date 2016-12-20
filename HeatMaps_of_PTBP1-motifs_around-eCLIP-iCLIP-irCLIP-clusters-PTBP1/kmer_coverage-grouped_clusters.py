'''
Created on Apr 25, 2014
 
@author: Nejc Haberman
 
description:
Script will accept fasta sequnces sorted by length. You need to specify a binning number maximum number of sequence length
where we look for the coverage based on all motifs (each motif must be in a new line).
For each sequnce we will look for the coverage (each nt overlapping with motif will count as 1) 
and sum together coverga from all fasta sequences.
 
input:
- input_fasta
- motifs (kmers)
- output_fasta
 
output:
- fasta with extra line of coverage (nt) per sequence
- total coverage 1 nt resolution of all fasta files 
'''
import sys
 
# load motifs from the file
def load_motifs(fin_fname_motifs):
    fin = open(fin_fname_motifs, "rt")
    motifs = []
    line = fin.readline()
    while line:
        motif = line.rstrip('\n').rsplit('\t')
        motif = str(motif[0]).upper()
        motifs.append(motif)
        line = fin.readline()
    fin.close()
    return motifs
 
# set the coverage of the sequence
def set_coverage(seq, motif, coverage):
    length = motif.__len__()
    motif_pos = seq.find(motif)
    pos = motif_pos
    while motif_pos != -1:
        coverage[pos:pos+length] = length * [1]
        motif_pos = seq[pos+1:].find(motif) #we search for the same motif downstream fro mthe previous one
        pos = motif_pos + pos + 1
    return coverage

def sum_coverage(total_coverage, adding_coverage):
    for i in range(0, total_coverage.__len__()):
        total_coverage[i] = total_coverage[i] + adding_coverage[i]
    return total_coverage 

# each number in array will be divided by a number
def coverage_norm(coverage, number):
    for i in range(0, coverage.__len__()):
        coverage[i] = float(coverage[i]) / float(number)
    return coverage 
 
def get_coverage(fin_fname, fin_fname_motifs, fout_fname):
    binning = 300
    max_length = 250
    fin = open(fin_fname, "rt")
    fout = open(fout_fname, "w")
    motifs = load_motifs(fin_fname_motifs)
    coverage = None
    total_coverage = None
    seq_length = 0
    line = fin.readline()
    info = ""
    counter = 0
    min_seq = 0
    max_seq = 0
    while line and int(max_seq) <= max_length:
        if line[0] == '>':
            info = line.rstrip('\n')
            counter += 1
        else:
            seq = line.rstrip('\n')
            seq = str(seq).upper()
            if seq_length == 0: #first time we need to get the length of the sequence
                seq_length = seq.__len__()
                min_seq = info.replace('>','')
                max_seq = info.replace('>','')
                coverage = [0] * max_length
                total_coverage = [0] * max_length
                
            if counter >= binning:
                max_seq = info.replace('>','')
                total_coverage = coverage_norm(total_coverage, counter-1) #divide array by counts (number of sequences) -1 because there is already a new one
                total_coverage[50-1] = 0.0000
                end = int(min_seq) + int((int(max_seq) - int(min_seq)) / 2.0)   #find end of the cluser position (if cluster is between 100..120, end will be 110)
                total_coverage[50+end+1] = 0.0000
                fout.write(str(min_seq) + '..' + str(max_seq) + "nt," + str(total_coverage).replace('[','').replace(']','') + '\n')
                min_seq = info.replace('>','')
                counter = 0
                total_coverage = [0] * max_length
             
            for i in range(0,motifs.__len__()): #for each motif we do a search
                motif = motifs[i]
                coverage = set_coverage(seq, motif, coverage)
            total_coverage = sum_coverage(total_coverage, coverage)
            coverage = [0] * max_length #initialize
        line = fin.readline()
    
    max_seq = info.replace('>','')
    total_coverage[50-1] = 0.0000
    #end = int(min_seq) + int((int(max_seq) - int(min_seq)) / 2.0)   #find end of the cluser position (if cluster is between 100..120, end will be 110)
    #total_coverage[50+end] = 0.0000
    total_coverage = coverage_norm(total_coverage, counter) #divide array by counts (number of sequences)
    fout.write(str(min_seq) + '..' + str(max_seq) + "nt," + str(total_coverage).replace('[','').replace(']','') + '\n')
    fin.close()
    fout.close()
               
'''
fin_fname_fasta = "/media/skgthab/SAMSUNG/UCL-backup/2014.06.09@Kulozik-short-long-iCLIPalignment/PTB/4.pipeline/heatmap_extended_clusters/filter2/Clusters-PTB-iCLIP-Er-merged-filtered50p-length-sorted-200nt-flanked_5_200.fasta"
fin_fname_motifs = "/media/skgthab/SAMSUNG/UCL-backup/2014.06.09@Kulozik-short-long-iCLIPalignment/PTB/4.pipeline/heatmap_extended_clusters/filter2/k-mers.tab"
fout_fname_fasta = "/media/skgthab/SAMSUNG/UCL-backup/2014.06.09@Kulozik-short-long-iCLIPalignment/PTB/4.pipeline/heatmap_extended_clusters/filter2/Clusters-PTB-iCLIP-Er-merged-filtered50p-length-sorted-200nt-flanked_5_200-HeatMap.csv"
get_coverage(fin_fname_fasta, fin_fname_motifs,fout_fname_fasta)
 
'''
 
if sys.argv.__len__() == 4:
    fin_fname_fasta = sys.argv[1]
    fin_fname_motifs = sys.argv[2]
    fout_fname_fasta = sys.argv[3]
    get_coverage(fin_fname_fasta, fin_fname_motifs, fout_fname_fasta)
else:
    #print str(sys.argv.__len__())
    print "error:\t2 arguments are needed\n" + '\n' +"example:\t $ python k-mer_coverage.py input_fname.fasta motifs.tab"
    
    
