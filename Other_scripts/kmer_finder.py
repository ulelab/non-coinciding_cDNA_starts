'''
Created on Aug 25, 2016

@author: Nejc Haberman

Input:    Fasta file
Output:   Text report of kmers frequency

Set max_kmer and min_kmer sequence that you are interested. By default it is set to 4 to search for tetramer sequences.

'''

import sys
import operator

# method will add kmer sequence to a dictionary or added counting +1 to existing one
def add_kmer(kmers, seq):
    if kmers.get(seq) == None:
        kmers.setdefault(seq, 1)
    else:
        kmers[seq] += 1
    return kmers

def find_kmers(fin_fname, fout_fname):
    fin = open(fin_fname, "rt")
    fout = open(fout_fname, "w")
    fout_end = open(fout_fname, "w")
    line = fin.readline()
    kmers = {}
    counter = 0 #number of sequences
    max_kmer = 4    #maximum kmer length
    min_kmer = 4    #minimum kmer length
    while line:
        if line[0] != '>':
            seq = line.rstrip('\n').upper()
            length = seq.__len__()
            for i in range(min_kmer, max_kmer+1):   #for each size of kmers
                start = 0
                while (start+i) <= length:  #we move through the sequence and add kmers to a dictionary
                    kmer = seq[start:start+i]
                    kmers = add_kmer(kmers, kmer)
                    start += 1
            counter += 1 
        line = fin.readline()
    sorted_kmers = sorted(kmers.items(), key=operator.itemgetter(1))    #sorting dictionary into array
    
    fout.write("kmer\tcount\tpercentage_of_inclusion\n")
    for i in range(sorted_kmers.__len__()-1, -1, -1):
        kmer, freq = sorted_kmers[i]
        perc = (float(freq) / float(counter)) * 100.0
        fout.write(kmer + "\t" + str(freq) + '\t' + str(perc) + '\n')
    fout.close()
    fout.close()
   
# main
if sys.argv.__len__() == 3:
    fin_fname = sys.argv[1]
    fout_fname = sys.argv[2]
    find_kmers(fin_fname, fout_fname)
else:
    #print str(sys.argv.__len__())
    print "error:\t2 arguments are needed\n" 

