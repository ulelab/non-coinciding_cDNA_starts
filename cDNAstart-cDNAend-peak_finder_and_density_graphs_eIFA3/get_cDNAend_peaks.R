'''
 This script returns peak positions from cDNA ends
 Input data are filtered cDNAs mapped to the transcript 
  args[1] <- "path to mapped cDNAs"
 Output data are selecred peaks
  args[2] <- "path to final cDNA end peak table"
'''

args<-commandArgs(TRUE)
reads <- read.table(args[1], sep='\t')
reads$V4 <- 1

# get star/end read positions related to exon exon junction
attach(reads)
reads$start <- V8 - V2  #all of them are on the plus strand and V9 you can use for read-end RNAmap
reads$end <- V9 - V2
detach(reads)

# get map of read end possitions for each exons
maps.end <- aggregate(V4 ~ end+V7+V2+V3, data=reads, FUN=sum) #we sum together normalized cDNA counts for each position
maps.end$junction <- paste(maps.end$V7, maps.end$V2, maps.end$V3, sep=':')
colnames(maps.end)[5] <- "end.peak"

maps.end <- maps.end[which(maps.end$end >= -25 & maps.end$end <= 20),]  #select end peaks in region -20 +20 from junction position
maps.end.max <- aggregate(end.peak ~ junction, data=maps.end, FUN=max)  #find a maximum peak
maps.selected <- merge(maps.end, maps.end.max, by=c("junction","end.peak"), all.y=TRUE) #select only maximum peaks

maps.selected2 <- maps.selected[which(maps.selected$end.peak >= median(maps.selected$end.peak)),] #filter peaks that have less then median cDNAs
maps.selected.avg <- aggregate(end ~ junction+end.peak, data=maps.selected2, FUN=mean) #in case of multiple max peaks select average position between them
maps.selected.avg$end <- as.integer(maps.selected.avg$end)  #convert average position to integer

# save the final sub set in tab delimited format
maps.selected.avg$junction.id <- lapply(strsplit(as.character(maps.selected.avg$junction), "\\:"), "[", 1)
maps.selected.avg$junction.start <- lapply(strsplit(as.character(maps.selected.avg$junction), "\\:"), "[", 2)
maps.selected.avg$junction.end <- lapply(strsplit(as.character(maps.selected.avg$junction), "\\:"), "[", 3)
final <- maps.selected.avg[c("junction.id","junction.start","junction.end","end.peak")]
write.table(final, args[2], row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
