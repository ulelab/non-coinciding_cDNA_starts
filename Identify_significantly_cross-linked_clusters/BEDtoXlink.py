'''
Created on Sep 16, 2016

@author: Nejc Haberman


Script will convert BED format read (after bedtools bamtobed) to cross link format (the beginning of the read)
'''

import sys

def BEDtoXlink(fname_in, fname_out):
    fin = open(fname_in, "rt")
    fout = open(fname_out, "w")
    line = fin.readline()

    while line:
        col = line.rstrip('\n').rsplit('\t')
        chr = col[0]
        start = int(col[1])
        end = int(col[2])
        strand = col[5]
        if strand == '+':
            fout.write(chr + '\t' + str(start - 1) + '\t' + str(start) + '\t' + '' + '\t' + '' + '\t' + strand + '\n')
        elif strand == '-':
            fout.write(chr + '\t' + str(end) + '\t' + str(end+1) + '\t' + '' + '\t' + '' + '\t' + strand + '\n')
        line = fin.readline()
    fout.close()
    fin.close()


if sys.argv.__len__() == 3:
    fname_in = sys.argv[1]
    fname_out = sys.argv[2]
    BEDtoXlink(fname_in, fname_out)
else:
    print("python BEDtoXlink.py <input_file> <output_file>")

