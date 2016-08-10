'''
Created on Aug 20, 2013

@author: Nejc

The script will flank the region in both directions in a new bed file.
'''

import sys

def flank_positions(fin_fname, fout_fname, left_shift, right_shift):
    fin = open(fin_fname, "rt")
    fout = open(fout_fname, "w")
    line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
        if line.__len__() > 6:
            chr = col[0]
            start = col[1]
            end = col[2]
            strand = col[5]
            if strand == '+':
                shift_start = int(start) - int(left_shift)
                shift_end = int(end) + int(right_shift)
            else:
                shift_end = int(end) + int(left_shift)
                shift_start = int(start) - int(right_shift)
            if shift_start >= 0 and shift_start < shift_end:
                fout.write(chr + '\t' + str(shift_start) + '\t' + str(shift_end) + '\t' + col[3] + '\t' + col[0] + ':' + col[1] + ':' + col[2] + '\t' + col[5] + '\n')
        line = fin.readline()

if sys.argv.__len__() == 5:
    fin_fname = sys.argv[1]
    fout_fname = sys.argv[2]
    left_shift = int(sys.argv[3])
    right_shift = int(sys.argv[4])
    flank_positions(fin_fname, fout_fname, left_shift, right_shift)
else:
    #print str(sys.argv.__len__())
    print "error:\t4 arguments are needed\n" + '\n' +"example:\t $ python bed_expand_positions.py input_fname.bed output_fname.bed left_shiftNUM right_shiftNUM"
