'''
Created on Sep 16, 2016

@author: Nejc Haberman


Script will sum together overlapping BED and added a count to a 5th column. BED needs to be sorted by gneomic positions and strand: sort -k1,1 -k2,2n -k6,6
'''

import sys

def BEDsum(fname_in, fname_out):
    fin = open(fname_in, "rt")
    fout = open(fname_out, "w")
    line = fin.readline()
    lastChr = None
    lastStart = None
    lastEnd = None
    lastStrand = None
    cDNAsum = 0
    while line:
        col = line.rstrip('\n').rsplit('\t')
        chr = col[0]
        start = col[1]
        end = col[2]
        cDNA = 1
        strand = col[5]
        if lastChr == None or (lastChr == chr and lastStart == start and lastEnd == end and lastStrand == strand):
            cDNAsum += int(cDNA)
        else:
            fout.write(lastChr + '\t' + lastStart + '\t' +  lastEnd + '\t' + '.' + '\t' + str(cDNAsum) + '\t' + lastStrand + '\n')
            cDNAsum = int(cDNA)
        lastChr = chr
        lastStart = start
        lastEnd = end
        lastStrand = strand
        line = fin.readline()
    fout.write(lastChr + '\t' + lastStart + '\t' +  lastEnd + '\t' + '.' + '\t' + str(cDNAsum) + '\t' + lastStrand + '\n')
    fout.close()
    fin.close()


if sys.argv.__len__() == 3:
    fname_in = sys.argv[1]
    fname_out = sys.argv[2]
    BEDsum(fname_in, fname_out)
else:
    print("python BEDsum.py <input_file> <output_file>")

