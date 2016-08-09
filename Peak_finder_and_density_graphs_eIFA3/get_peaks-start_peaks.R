'''
 This script returns peak positions from cDNA starts
 Input data are filtered cDNAs mapped to the transcript 
  args[1] <- "path to mapped cDNAs"
 Output data are selecred peaks
  args[2] <- "path to final cDNA start peak table"
'''

args<-commandArgs(TRUE)
reads <- read.table(args[1], sep='\t')
reads$V4 <- 1 #each cDNA count as one

# get star/end read positions related to exon exon junction
attach(reads)
reads$start <- V8 - V2  #all of them are on the plus strand and V9 you can use for read-end RNAmap
reads$end <- V9 - V2
detach(reads)

# get map of read start possitions for each exons
maps.start <- aggregate(V4 ~ start+V7+V2+V3, data=reads, FUN=sum) #we sum together normalized cDNA counts for each position
maps.start$junction <- paste(maps.start$V7, maps.start$V2, maps.start$V3, sep=':')
colnames(maps.start)[5] <- "start.peak"

maps.start <- maps.start[which(maps.start$start >= -100 & maps.start$start <= 0),]  #select end peaks in region -20 +20 from junction position
maps.start.max <- aggregate(start.peak ~ junction, data=maps.start, FUN=max)  #find a maximum peak
maps.selected <- merge(maps.start, maps.start.max, by=c("junction","start.peak"), all.y=TRUE) #select only maximum peaks

maps.selected2 <- maps.selected[which(maps.selected$start.peak >= median(maps.selected$start.peak)),] #filter peaks that have less then median cDNAs
maps.selected.avg <- aggregate(start ~ junction+start.peak, data=maps.selected2, FUN=mean) #in case of multiple max peaks select average position between them
maps.selected.avg$start <- as.integer(maps.selected.avg$start)  #convert average position to integer

# save the final sub set in tab delimited format
maps.selected.avg$junction.id <- lapply(strsplit(as.character(maps.selected.avg$junction), "\\:"), "[", 1)
maps.selected.avg$junction.start <- lapply(strsplit(as.character(maps.selected.avg$junction), "\\:"), "[", 2)
maps.selected.avg$junction.end <- lapply(strsplit(as.character(maps.selected.avg$junction), "\\:"), "[", 3)
final <- maps.selected.avg[c("junction.id","junction.start","junction.end","start.peak")]
write.table(maps.selected.avg, args[2], row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
