'''
Created on Oct 30, 2014

@author: Nejc Haberman

Mapped BED file will be converted into RNA map format with exon end position and crosslink around 100nt surrounding region.


'''

import sys

#script will filter exons that are less then 100 nt appart
def filter_exons(exons): 
    exons = map(int, exons)
    exons.sort()
    exons = exons[1:-2]     #last position of the exon end is ignored because there is no EJC complex
    filtered = []
    #print str(exons)
    filtered.append(exons[0])   #we always add the first one (there is no upstream exon position)
    for i in range(1, exons.__len__()):
        if abs(int(exons[i-1]) - int(exons[i])) < 100:  #ignore exon endings that have <100nt upstream exon ending
            #print str(exons[i]) + '\t' + str(exons[i-1])
            exons[i] = exons[i-1]   #previus exons stays the same because current one will be ignored
            #print str(exons)
        else:
            filtered.append(exons[i])
    return filtered 
    

def convert (fname_in, fname_out):
    fin = open(fname_in, "rt")
    fout = open(fname_out, "w")
    line = fin.readline()
    while line:
        #clean empty spaces
        line = line.replace("           ","\t").replace("          ","\t").replace("         ","\t").replace("        ","\t").replace("       ","\t").replace("      ","\t").replace("     ","\t").replace("    ","\t").replace("    ","\t").replace("   ","\t").replace("  ","\t").replace(" ","\t")
        col = line.rstrip('\n').rsplit("\t")
        #print line
        strand = col[-1]
        if strand == '+':   #minus strand reads are ignored because reads were mapped to a single transcript
            gene = col[0]
            gene_info = gene.rsplit('|')
            gene_id = gene_info[0]
            #print line
            exon_ends = gene_info[2].rsplit(';')
            if exon_ends.__len__() > 4: #we ignore transcripts with less then 4 exons ()
                exon_ends = filter_exons(exon_ends)
                read_start = int(col[1])
                read_end = int(col[2])
                
                for i in range(0, exon_ends.__len__()): #we save read position to each exonic endings
                    ref_pos = exon_ends[i]
                    if (abs(int(ref_pos) - read_start) <= 250): #we want draw a map that is more then 250 nt appart 
                        fout.write(gene_id + '\t' + str(ref_pos) + '\t' + str(int(ref_pos)+1) + '\t\t\t' + '+' + '\t' + gene_id + '\t' + str(read_start) + '\t' + str(read_end) + '\t\t\t' + '+' + '\n')
        line = fin.readline()
    fout.close()
    fin.close()


'''
convert("/media/skgthab/SAMSUNG/UCL-backup/2015.04.10@eIF4E3-Zhen/mapping_to_transcripts/test/2012_0003_NoIndex_L002_R1_001._CAAT.fq-2mis-uniq.bed", "/media/skgthab/SAMSUNG/UCL-backup/2015.04.10@eIF4E3-Zhen/mapping_to_transcripts/test/2012_0003_NoIndex_L002_R1_001._CAAT.fq-2mis-uniq-RNAmap.bed")
'''
    
if sys.argv.__len__() == 3:
    fname_in = sys.argv[1]
    fname_out = sys.argv[2]
    convert(fname_in, fname_out)
else:
    print("python toOneLine.py <input_file> <output_file>")
