'''
 This script filters out only top 1000 transcripts by the number of cDNAs
 Input data is from mapped cDNAs to the transcriptome sequence for eIFA3 RNA Binding Protein:
  args[1] <- "path to mapped cDNAs"
 Output data are filtered cDNAs from top 1000 transcripts
  args[2] <- "path to filtered cDNAs"
'''

args<-commandArgs(TRUE)
reads <- read.table(args[1], sep='\t')
reads$V4 <- 1 #each cDNA count as one
read.number <- nrow(reads)  #get number of all cDNAs

#select top 1000 transcripts
transcripts <- aggregate(V4 ~ V1, data=reads, FUN = sum)  #get a number of cDNAs for each transcript
colnames(transcripts)[2] <- "transcript.cDNA.count"
transcripts <- transcripts[order(-transcripts$transcript.cDNA.count),]
top1000.transcripts <- transcripts[0:1000,]
reads <- merge(reads, top1000.transcripts, by="V1")
reads$V4 <- 1 / read.number #normalize

# see how many normalized reads per exon do we have
exons <- aggregate(V4 ~ V13, data=reads, FUN=sum)
exons <- exons[order(exons$V4),]

# now we select reads from selected exons
selected.reads <- merge(reads, exons, by="V13", all.y=TRUE)
selected.reads$V1 <- NULL
write.table(selected.reads, args[2], row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
