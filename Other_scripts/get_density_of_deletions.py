'''
Created on June 25, 2016
 
@author: Nejc Haberman
 
description:
Script will recive a SAM file format and return a density array of positions where deletions are across cDNAs
For more details about SAM format: https://samtools.github.io/hts-specs/SAMv1.pdf
'''
import sys

def get_deletions_coverage (fin_sam):
    fin = open(fin_sam, "rt")
    line = fin.readline()
    deletion_sum = 0    #number of all deletions
    read_count = 0
    coverage = 150 * [0] #max cDNA length set to 150 nt
    deletion_reads_count = 0
    while line[0] == '@':    #ignore header of .SAM
        line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
        strand = col[1] # 0 is plus strand, 16 is minus
        cigar = col[5]
        seq = col[9]
        length = seq.__len__()
        deletion_pos = cigar.find("D")  #D flag in cigar SAM format stands for a deletion 
        read_count += 1
        if deletion_pos != -1:
            deletion_reads_count += 1
            deletions = cigar.rsplit("D")   #thre can be multiple deletions per cDNA
            for i in range(0, cigar.count("D")):    #loop is based on the number of deletions
                deletion = deletions[i]
                deletion_num = int(deletion[-1])    #number of deletion is always at last char before D which we used for spliting
                deletion_sum += deletion_num
                upstream_positions = deletion[0:-1].replace("M",",").replace("I",",").replace("'","").rstrip(",")
                upstream_positions = upstream_positions.rsplit(",")
                upstream_positions_array = map(int, upstream_positions)
                upstream_distance = sum(upstream_positions_array)
                if strand == "0":   #plus strand
                    for i in range(0,deletion_num): #if we have a multiple deletions they can also be one after another
                        coverage[upstream_distance+i] += 1
                        if upstream_distance <= 8 and deletion_num == 1:
                            fout_up.write(line)
                        else:
                            fout_down.write(line)
                if strand == "16":  #minus strand
                    coverage = coverage[::-1]   #we reverse it
                    for i in range(0,deletion_num):
                        coverage[upstream_distance+i] += 1
                        if (length - upstream_distance) <= 8 and deletion_num == 1:
                            fout_up.write(line)
                        else:
                            fout_down.write(line)
                    coverage = coverage[::-1]   #and reverse it to original 
        line = fin.readline()
        
    print "\nNumber of all deletions: " + str(deletion_sum) + " (" + str(float(deletion_sum)/float(read_count) * 100.0) + " %)" 
    print "Number of cDNAs with deletions: " + str(deletion_reads_count) + " (" + str(float(deletion_reads_count)/float(read_count) * 100) + " %)" 
    print "Number of all cDNAs: " + str(read_count) + '\n'
    print "\nPositions of all deletions across cDNAs:\n"
    print str(coverage).replace('[','').replace(']','').replace("'","")
    fin.close()

# main
if sys.argv.__len__() == 2:
    fin_sam = sys.argv[1]
    get_deletions_coverage(fin_sam)
else:
    print "you need 1 argument to run the script in SAM file format"
    quit()

    
    