'''
Created on Nov 28, 2013

@author: Nejc Haberman


The script will read fasta file and remove random barcode and experimental barcode from fasta. Random barcode will be saved in the header of fasta file
'''


import sys

def swap_barcodes(fin_fasta, fout_fasta):
    finFasta = open(fin_fasta, "rt")
    foutFasta = open(fout_fasta, "w")
    line = finFasta.readline()
    while line:
        if line[0] != '>':
            randomBarcode = line[0:3] + line[7:9]
            experimentBarcode = line[3:7]
            seqRead = line[9:]

            foutFasta.write(">"+randomBarcode + '\n')
            foutFasta.write(seqRead)
        line = finFasta.readline()
    finFasta.close()
    foutFasta.close()

if sys.argv.__len__() == 3:
    fin_fasta = sys.argv[1]
    fout_fasta = sys.argv[2]
    swap_barcodes(fin_fasta, fout_fasta)
else:
    print "you need 2 arguments to run the script"
    quit()

