library("ggplot2")
library("smoother")

args<-commandArgs(TRUE)
smoothing_window <- 5

#colour-9
colour30 <- "#242B38"
colour35 <- "#2165E8"
colour40 <- "#2CA3B5"
colour60 <- "#EA882C"
colour130 <- "#BF1F1B"
colourLong <- "#1D912F"  #untrimmed

####################

# TO DO
# add table of selected junctions and merge them with the ones here V2 and V3

#RNAmap.1 <- read.table("eIF4E3_4A3_all-complete-filter100.bed", sep='\t')
#RNAmap.1.long <- read.table("/media/skgthab/storage_backup/UCL/2015.08.05@eIF4A3-Ina/exon-classification-Ina-2/RNAmaps/iCLIP_EIF4A3_293Flp_RNase1:500_Hs_NNNGAGCNN_20150803_Mx1_17_18.fq-truncated.fq-incomplete-unmapped.fasta.trimmed-2mis-uniq-merged-filtered100.bed", sep='\t')
#junctions <- read.table("eIF4E3_4A3_all-complete.bed-filtered100.bed-filtered100-selected-small-selected-max_peaks.tab", sep='\t')
junctions <- read.table(args[3])

RNAmap.1 <- read.table(args[1], sep='\t')

# select only reads that are from selected junctions
RNAmap.1 <- merge(RNAmap.1, junctions, by=c("V1","V2","V3"),all=FALSE)
RNAmap.1$V2 <- RNAmap.1$V2 + RNAmap.1$V4.y  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! here we change junction position to max peak position !!!!!!!!!!!!!!!!!!!!!!!!!!
RNAmap.1$V4 <- 1  #count
RNAmap.1$V10 <- RNAmap.1$V9 - RNAmap.1$V8

RNAmap.1.long <- read.table(args[2], sep='\t')
RNAmap.1.long <- merge(RNAmap.1.long, junctions, by=c("V1","V2","V3"),all=FALSE) 
RNAmap.1.long$V2 <- RNAmap.1.long$V2 + RNAmap.1.long$V4.y  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! here we change junction position to max peak position !!!!!!!!!!!!!!!!!!!!!!!!!!
RNAmap.1.long$V4 <- 1  #count
RNAmap.1.long$V10 <- RNAmap.1.long$V9 - RNAmap.1.long$V8

RNAmap.1$start <- RNAmap.1$V8 - RNAmap.1$V2  #all of them are on the plus strand and V9 you can use for read-end RNAmap
RNAmap.1$end <- RNAmap.1$V9 - RNAmap.1$V2

RNAmap.1.long$start <- RNAmap.1.long$V8 - RNAmap.1.long$V2  #all of them are on the plus strand and V9 you can use for read-end RNAmap
RNAmap.1.long$end <- RNAmap.1.long$V9 - RNAmap.1.long$V2

#normalize by the number of reads and exon junctions 
RNAmap.1$norm <- RNAmap.1$V4 / nrow(RNAmap.1[!duplicated(RNAmap.1[2:3]),])
RNAmap.1.long$norm <- RNAmap.1.long$V4 / nrow(RNAmap.1.long[!duplicated(RNAmap.1.long[2:3]),])

#normalize by the total mRNA coverage
mRNAstarts <- aggregate(V2 ~ V1, data=RNAmap.1, FUN=min)
mRNAends <- aggregate(V2 ~ V1, data=RNAmap.1, FUN=max)
mRNA <- merge(mRNAstarts, mRNAends, by="V1")
mRNA$length <- mRNA$V2.y - mRNA$V2.x
RNAmap.1$norm <- RNAmap.1$norm * sum(mRNA$length)

mRNAstartsLong <- aggregate(V2 ~ V1, data=RNAmap.1.long, FUN=min)
mRNAendsLong <- aggregate(V2 ~ V1, data=RNAmap.1.long, FUN=max)
mRNAlong <- merge(mRNAstartsLong, mRNAendsLong, by="V1")
mRNAlong$length <- mRNAlong$V2.y - mRNAlong$V2.x
RNAmap.1.long$norm <- RNAmap.1.long$norm * sum(mRNAlong$length)


RNAmap.1.map.start30 <- RNAmap.1[which(RNAmap.1$V10 < 30),]
RNAmap.1.map.start30$norm <- RNAmap.1.map.start30$norm / sum(RNAmap.1.map.start30$V4)
RNAmap.1.map.end30 <- RNAmap.1[which(RNAmap.1$V10 < 30),]
RNAmap.1.map.end30$norm <- RNAmap.1.map.end30$norm / sum(RNAmap.1.map.end30$V4)
RNAmap.1.map.start30 <- aggregate(norm ~ start, data=RNAmap.1.map.start30, FUN=sum) #we sum together normalized cDNA counts for each position
RNAmap.1.map.end30 <- aggregate(norm ~ end, data=RNAmap.1.map.end30, FUN=sum) #we sum together normalized cDNA counts for each position
colnames(RNAmap.1.map.start30)[2] <- "start30"
colnames(RNAmap.1.map.end30)[2] <- "end30"

RNAmap.1.map.start35 <- RNAmap.1[which(RNAmap.1$V10 >= 30 & RNAmap.1$V10 < 35),]
RNAmap.1.map.start35$norm <- RNAmap.1.map.start35$norm / sum(RNAmap.1.map.start35$V4)
RNAmap.1.map.end35 <- RNAmap.1[which(RNAmap.1$V10 >= 30 & RNAmap.1$V10 < 35),]
RNAmap.1.map.end35$norm <- RNAmap.1.map.end35$norm / sum(RNAmap.1.map.end35$V4)
RNAmap.1.map.start35 <- aggregate(norm ~ start, data=RNAmap.1.map.start35, FUN=sum) #we sum together normalized cDNA counts for each position
RNAmap.1.map.end35 <- aggregate(norm ~ end, data=RNAmap.1.map.end35, FUN=sum) #we sum together normalized cDNA counts for each position### set all values to 0 in range -500.500
colnames(RNAmap.1.map.start35)[2] <- "start35"
colnames(RNAmap.1.map.end35)[2] <- "end35"

RNAmap.1.map.start40 <- RNAmap.1[which(RNAmap.1$V10 >= 35 & RNAmap.1$V10 < 40),]
RNAmap.1.map.start40$norm <- RNAmap.1.map.start40$norm / sum(RNAmap.1.map.start40$V4)
RNAmap.1.map.end40 <- RNAmap.1[which(RNAmap.1$V10 >= 35 & RNAmap.1$V10 < 40),]
RNAmap.1.map.end40$norm <- RNAmap.1.map.end40$norm / sum(RNAmap.1.map.end40$V4)
RNAmap.1.map.start40 <- aggregate(norm ~ start, data=RNAmap.1.map.start40, FUN=sum) #we sum together normalized cDNA counts for each position
RNAmap.1.map.end40 <- aggregate(norm ~ end, data=RNAmap.1.map.end40, FUN=sum) #we sum together normalized cDNA counts for each position### set all values to 0 in range -500.500
colnames(RNAmap.1.map.start40)[2] <- "start40"
colnames(RNAmap.1.map.end40)[2] <- "end40"

RNAmap.1.map.start60 <- RNAmap.1[which(RNAmap.1$V10 >= 40),]
RNAmap.1.map.start60$norm <- RNAmap.1.map.start60$norm / sum(RNAmap.1.map.start60$V4)
RNAmap.1.map.end60 <- RNAmap.1[which(RNAmap.1$V10 >= 40),]
RNAmap.1.map.end60$norm <- RNAmap.1.map.end60$norm / sum(RNAmap.1.map.end60$V4)
RNAmap.1.map.start60 <- aggregate(norm ~ start, data=RNAmap.1.map.start60, FUN=sum) #we sum together normalized cDNA counts for each position
RNAmap.1.map.end60 <- aggregate(norm ~ end, data=RNAmap.1.map.end60, FUN=sum) #we sum together normalized cDNA counts for each position### set all values to 0 in range -500.500
colnames(RNAmap.1.map.start60)[2] <- "start60"
colnames(RNAmap.1.map.end60)[2] <- "end60"

RNAmap.1.long$norm <- RNAmap.1.long$norm / sum(RNAmap.1.long$V4)
RNAmap.1.map.startLong <- aggregate(norm ~ start, data=RNAmap.1.long, FUN=sum) #we sum together normalized cDNA counts for each position
RNAmap.1.map.endLong <- aggregate(norm ~ end, data=RNAmap.1.long, FUN=sum) #we sum together normalized cDNA counts for each position### set all values to 0 in range -500.500
colnames(RNAmap.1.map.startLong)[2] <- "startLong"
colnames(RNAmap.1.map.endLong)[2] <- "endLong"

zero <- 0
zero$start <- c(-250:250)
RNAmap.1.map.start <- merge(zero, RNAmap.1.map.start30, by="start",all=TRUE)
RNAmap.1.map.start <- merge(RNAmap.1.map.start, RNAmap.1.map.start35, by="start",all=TRUE)
RNAmap.1.map.start <- merge(RNAmap.1.map.start, RNAmap.1.map.start40, by="start",all=TRUE)
RNAmap.1.map.start <- merge(RNAmap.1.map.start, RNAmap.1.map.start60, by="start",all=TRUE)
RNAmap.1.map.start <- merge(RNAmap.1.map.start, RNAmap.1.map.startLong, by="start",all=TRUE)
RNAmap.1.map.start[is.na(RNAmap.1.map.start)]<-0

zero <- 0
zero$end <- c(-250:250)
RNAmap.1.map.end <- merge(zero, RNAmap.1.map.end30, by="end",all=TRUE)
RNAmap.1.map.end <- merge(RNAmap.1.map.end, RNAmap.1.map.end35, by="end",all=TRUE)
RNAmap.1.map.end <- merge(RNAmap.1.map.end, RNAmap.1.map.end40, by="end",all=TRUE)
RNAmap.1.map.end <- merge(RNAmap.1.map.end, RNAmap.1.map.end60, by="end",all=TRUE)
RNAmap.1.map.end <- merge(RNAmap.1.map.end, RNAmap.1.map.endLong, by="end",all=TRUE)
RNAmap.1.map.end[is.na(RNAmap.1.map.end)]<-0

# smoothing
RNAmap.1.map.start$smooth30 <- smth(RNAmap.1.map.start$start30, window = smoothing_window, method = "gaussian") #SMOOTHING
RNAmap.1.map.start$smooth35 <- smth(RNAmap.1.map.start$start35, window = smoothing_window, method = "gaussian") #SMOOTHING
RNAmap.1.map.start$smooth40 <- smth(RNAmap.1.map.start$start40, window = smoothing_window, method = "gaussian") #SMOOTHING
RNAmap.1.map.start$smooth60 <- smth(RNAmap.1.map.start$start60, window = smoothing_window, method = "gaussian") #SMOOTHING
RNAmap.1.map.start$smoothLong <- smth(RNAmap.1.map.start$startLong, window = smoothing_window, method = "gaussian") #SMOOTHING

RNAmap.1.map.end$smooth30 <- smth(RNAmap.1.map.end$end30, window = smoothing_window, method = "gaussian") #SMOOTHING
RNAmap.1.map.end$smooth35 <- smth(RNAmap.1.map.end$end35, window = smoothing_window, method = "gaussian") #SMOOTHING
RNAmap.1.map.end$smooth40 <- smth(RNAmap.1.map.end$end40, window = smoothing_window, method = "gaussian") #SMOOTHING
RNAmap.1.map.end$smooth60 <- smth(RNAmap.1.map.end$end60, window = smoothing_window, method = "gaussian") #SMOOTHING
RNAmap.1.map.end$smoothLong <- smth(RNAmap.1.map.end$endLong, window = smoothing_window, method = "gaussian") #SMOOTHING

tans <- 0.8
g_size <- 0.4
ggplot() + theme_bw() +
  geom_line(aes(RNAmap.1.map.start$start, as.vector(RNAmap.1.map.start$smooth60)),size=g_size, colour = colour60, alpha=tans) + 
  geom_line(aes(RNAmap.1.map.start$start, as.vector(RNAmap.1.map.start$smooth40)),size=g_size, colour = colour40, alpha=tans) + 
  geom_line(aes(RNAmap.1.map.start$start, as.vector(RNAmap.1.map.start$smooth35)),size=g_size, colour = colour35, alpha=tans) + 
  geom_line(aes(RNAmap.1.map.start$start, as.vector(RNAmap.1.map.start$smooth30)),size=g_size, colour = colour30, alpha=tans) + 
  geom_line(aes(RNAmap.1.map.end$end, as.vector(RNAmap.1.map.end$smooth60)),size=g_size, colour = colour60, alpha=tans, linetype=3) + 
  geom_line(aes(RNAmap.1.map.end$end, as.vector(RNAmap.1.map.end$smooth40)),size=g_size, colour = colour40, alpha=tans, linetype=3) + 
  geom_line(aes(RNAmap.1.map.end$end, as.vector(RNAmap.1.map.end$smooth35)),size=g_size, colour = colour35, alpha=tans, linetype=3) + 
  geom_line(aes(RNAmap.1.map.end$end, as.vector(RNAmap.1.map.end$smooth30)),size=g_size, colour = colour30, alpha=tans, linetype=3) + 
  ggtitle(paste(args[4],"- read starts", sep=" ")) + 
  xlab("position relative reference") + 
  ylab("normalized cDNA counts") + 
  theme(text=element_text(size=6),axis.text=element_text(size=6), axis.title=element_text(size=6,face="plain")) + 
  scale_x_continuous(limits = c(-75, 75)) #+ 
  #scale_y_continuous(limits = c(0, 3))



