# Non-coinciding cDNA starts

#Trimming of adapter sequences
Before mapping the cDNAs, we removed random barcodes and trimmed the 3´ Solexa adapter sequence. Adapter sequences were trimmed by FASTX-Toolkit 0.0.13 adapter removal software, using the following parameters: -Q 33 -a AGATCGGAAG -c -n -l 26. For reads that did not contain parts of the adapter sequence, -C parameter was used, and these were analyzed separately.

#Genome mapping of PTBP1 and U2AF65 iCLIP
We used UCSC hg19/GRCh37 genome assembly and bowtie2 2.1 alignment software with default settings accepting uniquely mapped cDNAs to a single genomic position and allowing maximum of 2 missmatches. After mapping, cDNAs with the same random barcode that mapped to the same starting position on the genome were considered to result from PCR amplification and were collapsed to a single cDNA.

#Transcriptome mapping of eIF4A3 iCLIP
For mapping, we compiled a set of representative mRNA sequences from BioMart Ensembl Genes 79, where we used the longest mRNA sequence available for each gene. We mapped to these mRNAs with Bowtie2.1 alignment software, allowing 2 mismatches. After mapping, cDNAs with the same random barcode that mapped to the same starting position on an mRNA were considered to result from PCR amplification and were collapsed to a single cDNA.

#Classification of cDNA length
Only cDNAs that mapped to a unique genomic position were evaluated. These were separated into cDNAs that were <30 nt, 30-34 nt, 35-39 nt or >40 long after trimming. 

#Definition of crosslink-associated motifs
We reasoned that sequence motifs enriched directly at the starts of the control eCLIP cDNAs might uncover preferences of UV crosslinking, since they are thought to represent a mixture of crosslink sites for many different RBPs, and thus they should not reflect sequence specificity of any specific RBP. We therefore examined occurrence of tetramers that overlapped with the nucleotide preceding the cDNA-starts in PTBP1 control iCLIP (position -1) in comparison with the ones overlapping with the 10th nucleotide preceding the cDNA-starts (position -10). Tetramers that are enriched over 1.5 fold at position -1 compared to -10 include. We excluded the TTTT tetramer from further analyses, since it is often part of longer tracts of Ts, and therefore its inclusion decreases the resolution of analysis. Thus, TTTG, TTTC, TTGG, TTTA, ATTG, ATTT, TCGT, TTGA, TTCT and CTTT were used for all analyses of crosslink-associated motifs.

#Definition of PTBP-target motifs
To identify the motifs bound by PTBP1, we searched for pentamers enriched in the region [-10..10] around the cDNA-start peaks identified in each crosslink cluster defined by PTBP1-iCLIP2.  69 pentamers had enrichment z-score > 299 and were used as PTBP1-target pentamers for further analyses. Their sequences are: TCTTT, CTTTC, TCTTC, CTTCT, TCTCT, CTCTC, TTTCT, TTCTC, TTCTT, TTTTC, TCCTT, CTCTT, ATTTC, TTCCT, CTTCC, TTTCC, CCTTT, CTTTT, CCTTC, TCTGT, TTCTG, TCCTC, CTTCA, ATCTT, TGTCT, TCTGC, CTCCT, CCTCT, GTCTT, TCTAT, TCTCC, ATTCC, TTCTA, CTTTG, TATCT, ACTTC, TTATC, CTTAT, CTATT, TTCAT, TTCCA, TCTTG, TTGTC, TTGCT, CTCTA, CTCTG, TATTT, TCCCT, TCATT, TTCCC, CATTT, ATTCT, TTTAC, GTTCT, CTATC, TCATC, CTTTA, TGTTC, TATTC, CATCT, TACTT, CTGTT, CTTGC, ACCTT, TTTCA, TTTGT, TGTTT, CTTGT, ACTTT. All of these pentamers are enriched in pyrimidines, in agreement with the known preference of PTBP1 for UC-rich binding motifs.

#Normalisation of data for drawing of density graphs
All normalisations were performed in R (version 3.1.0) together with “ggplot2” and “smoother” package for final graphical output.

- For analysis of EIFA3 iCLIP, each density graph (RNA map) shows a distribution of cDNA start and end positions relative to positions of exon-exon junctions in mRNAs. All exon-exon junctions within coding regions were taken into account, apart from those that junctioned to the first or last exon in the mRNA. The number of cDNAs starting or ending at each position on the graph was normalised by the number of all cDNAs mapped to representative mRNAs, the mRNA length and the number of examined exon-exon junction positions in the following way: 
RNAmap[n] = ((cDNAs[n] / sum(cDNAs)) * length(mRNAs) / count(exon_exon_junctions)
where [n] stands for a specific position on the density graph.
To draw the graph, we then used Gaussian method with a 3-nucleotide smoothing window. 

- For analysis of PTBP1, U2AF65 iCLIP and CLIP, each density graph (RNA map) shows a distribution of cDNA start and end positions relative to positions of its binding sites. We defined the binding sites in two ways: either by using the position of Y-tracts or extended crosslink clusters in introns. All introns in protein coding genes were taken into account. Due to highly variable abundance of intronic RNAs (and occasional presence of highly abundant non-coding transcripts, such as snoRNAs), we first divided counts at each binding site by dividing them by the count at the MaxCount. To define MaxCount, we examined the region of the binding site, as well as 120 nt 5´ and 3´ of the binding site, to find the nucleotide with the largest count of cDNA count starts or cDNA ends (according to whether starts or ends were plotted on the graph). This avoids from the highly abundant intronic snoRNAs or other abundant introns from dominating the results. The MaxCount-normalised counts of cDNAs starting or ending at each position on the graph was then normalised by the density of all MaxCount-normalised cDNA count/nt starting in the region 50-100 nt downstream of the binding site in the following way: 
RNAmap[n] = (MaxCount-normalised cDNAs[n]) / MaxCount-normalised cDNA density(outside the binding site)
where [n] stands for a specific position on the density graph.
To draw the graph, we then used Gaussian method with a 10-nucleotide smoothing window. 

#Assignment of the cDNA-end peak in eIF4A3 iCLIP
For cDNA-end peak assignment in eIFA3 iCLIP data, we used exons longer then 100 nt that were in the top 50% of the distribution of exons based on cDNA coverage. This ensured that sufficient cDNAs were available for assignment of the putative binding sites. We then summarised all cDNA-end positions in the region -20 to +25 around exon-exon junctions and selected the position with the maximum cDNA count as the ‘cDNA-end peak’.

#Analysis of pairing probability
Computational prediction of the secondary structure around the cDNA-end peaks was performed using the RNAfold program with the default parameters.

#Analysis of cDNA transitions
Density of U>T transitions across cDNAs was performed by using the samtools software with the the following parameters: samtools calmd –u –u genomic.fasta input_BAM > BAM-with_transitions. This pipeline replaces BAM format mapped cDNA sequences with transitions relative to genomic reference. This step was performed with a custom python script (available on github repository) that returns a density array of U>T transitions for cDNAs that are shorter then 40 nts. For the final visualisation of density graphs we used the same approach as for all other density figures. 

# Main scripts
 - mapping_to_genome-PTBP1_and_U2AF65_iCLIP-pipeline.sh (PTBP1 and U2AF65 mapping pipeline to genome)
 - mapping_to_genome-PTBP1_and_U2AF65_CLIP-pipeline.sh (PTBP1 and U2AF65 mapping pipeline to genome)
 - mapping_to_transcriptome-eIFA3_iCLIP-pipeline.sh (eIFA3 mapping pipeline to transcriptome)
 - mapping_to_transcriptome-eIFA3_CLIP-pipeline.sh (eIFA3 mapping pipeline to transcriptome)
 - normalisation_and_density_graph-eIFA4.R (normalisation and drowing of eIFA3)
 - normalisation_and_density_graph-PTBP1-U2AF65.sh (normalisation and drowing of PTBP1 and U2AF65)
 - main_get_cDNA_peaks.sh (eIFA3 selecting top 1000 transcripts and reporiting/drawing around cDNA start and cDNA end peaks)
 - other scripts (kmer finder, flanking BED positions, density of deletions across all cDNAs)
