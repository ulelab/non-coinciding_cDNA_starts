# Final result of the script is a normalised figure around exon ends showing a density of complete and incomplete reads.
# Input data are complete and incomplete reads mapped to the transcriptome sequence for eIFA3 RNA Binding Protein

library("ggplot2")
library("smoother")

args<-commandArgs(TRUE)
smoothing_window <- 3

# importing the data
RNAmap.complete <- read.table(args[1], sep='\t')
RNAmap.complete$V4 <- 1  #each read count as 1
RNAmap.complete$V10 <- RNAmap.complete$V9 - RNAmap.complete$V8  #read length

RNAmap.incomplete <- read.table(args[2], sep='\t')  #read length
RNAmap.incomplete$V4 <- 1  #each read count as 1
RNAmap.incomplete$V10 <- RNAmap.incomplete$V9 - RNAmap.incomplete$V8

# if you want to filter reads by read length
#RNAmap.complete <- RNAmap.complete[which(RNAmap.complete$V10 >=17 & RNAmap.complete$V10 < 30),]
#RNAmap.incomplete <- RNAmap.incomplete[which(RNAmap.incomplete$V10 >=17 & RNAmap.incomplete$V10 < 30),]

# setting positions
RNAmap.complete$start <- RNAmap.complete$V8 - RNAmap.complete$V2  #read start position relative to exon end
RNAmap.complete$end <- RNAmap.complete$V9 - RNAmap.complete$V2     #read end position relative to exon end
RNAmap.incomplete$start <- RNAmap.incomplete$V8 - RNAmap.incomplete$V2  
RNAmap.incomplete$end <- RNAmap.incomplete$V9 - RNAmap.incomplete$V2

#normalise by the number of reads and exon junctions 
RNAmap.complete$norm <- RNAmap.complete$V4 / sum(RNAmap.complete$V4)
RNAmap.incomplete$norm <- RNAmap.incomplete$V4 / sum(RNAmap.incomplete$V4)
RNAmap.complete$norm <- RNAmap.complete$norm / nrow(RNAmap.complete[!duplicated(RNAmap.complete[2:3]),])
RNAmap.incomplete$norm <- RNAmap.incomplete$norm / nrow(RNAmap.incomplete[!duplicated(RNAmap.incomplete[2:3]),])

#normalise by the total mRNA coverage
mRNAstarts <- aggregate(V2 ~ V1, data=RNAmap.complete, FUN=min)
mRNAends <- aggregate(V2 ~ V1, data=RNAmap.complete, FUN=max)
mRNA <- merge(mRNAstarts, mRNAends, by="V1")
mRNA$length <- mRNA$V2.y - mRNA$V2.x
RNAmap.complete$norm <- RNAmap.complete$norm * sum(mRNA$length)

mRNAstartsLong <- aggregate(V2 ~ V1, data=RNAmap.incomplete, FUN=min)
mRNAendsLong <- aggregate(V2 ~ V1, data=RNAmap.incomplete, FUN=max)
mRNAlong <- merge(mRNAstartsLong, mRNAendsLong, by="V1")
mRNAlong$length <- mRNAlong$V2.y - mRNAlong$V2.x
RNAmap.incomplete$norm <- RNAmap.incomplete$norm * sum(mRNAlong$length)

# sum cDNA counts for each position
RNAmap.complete.map.start <- aggregate(norm ~ start, data=RNAmap.complete, FUN=sum) #we sum together normalised cDNA counts for each read start position
RNAmap.complete.map.end <- aggregate(norm ~ end, data=RNAmap.complete, FUN=sum) #we sum together normalised cDNA counts for each read end position
RNAmap.incomplete.map.start <- aggregate(norm ~ start, data=RNAmap.incomplete, FUN=sum) #we sum together normalised cDNA counts for each read start position
RNAmap.incomplete.map.end <- aggregate(norm ~ end, data=RNAmap.incomplete, FUN=sum) #we sum together normalised cDNA counts for each read end position

# initialize and create a final map positions
zero <- 0
zero$start <- c(-250:250)
RNAmap.complete.map.start <- merge(zero, RNAmap.complete.map.start, by="start",all=TRUE)
RNAmap.complete.map.start[is.na(RNAmap.complete.map.start)]<-0
zero <- 0
zero$end <- c(-250:250)
RNAmap.complete.map.end <- merge(zero, RNAmap.complete.map.end, by="end",all=TRUE)
RNAmap.complete.map.end[is.na(RNAmap.complete.map.end)]<-0

zero <- 0
zero$start <- c(-250:250)
RNAmap.incomplete.map.start <- merge(zero, RNAmap.incomplete.map.start, by="start",all=TRUE)
RNAmap.incomplete.map.start[is.na(RNAmap.incomplete.map.start)]<-0
zero <- 0
zero$end <- c(-250:250)
RNAmap.incomplete.map.end <- merge(zero, RNAmap.incomplete.map.end, by="end",all=TRUE)
RNAmap.incomplete.map.end[is.na(RNAmap.incomplete.map.end)]<-0

# smoothing
RNAmap.complete.map.start$smooth <- smth(RNAmap.complete.map.start$norm, window = smoothing_window, method = "gaussian") 
RNAmap.complete.map.end$smooth <- smth(RNAmap.complete.map.end$norm, window = smoothing_window, method = "gaussian") 
RNAmap.incomplete.map.start$smooth <- smth(RNAmap.incomplete.map.start$norm, window = smoothing_window, method = "gaussian") 
RNAmap.incomplete.map.end$smooth <- smth(RNAmap.incomplete.map.end$norm, window = smoothing_window, method = "gaussian") 


tans <- 1
g_size <- 0.4
RNAmap.complete.map.gg.start <- ggplot() + theme_bw() +
  geom_line(aes(RNAmap.complete.map.start$start, as.vector(RNAmap.complete.map.start$smooth)),size=g_size, colour = "#242B38", cut=TRUE, fill = "#242B38", alpha=tans, adjust = 0.4) + 
  geom_line(aes(RNAmap.incomplete.map.start$start, as.vector(RNAmap.incomplete.map.start$smooth)),size=g_size, colour = "#1D912F", cut=TRUE, fill = "#1D912F", alpha=tans, adjust = 0.4) + 
  ggtitle("eIFA3 read start positions") + 
  xlab("position relative to exon end") + 
  ylab("normalised cDNA counts") + 
  theme(text=element_text(size=6),axis.text=element_text(size=6), axis.title=element_text(size=6,face="plain")) + 
  scale_x_continuous(limits = c(-75, 25))

RNAmap.complete.map.gg.end <- ggplot() + theme_bw() +
  geom_line(aes(RNAmap.complete.map.end$end, as.vector(RNAmap.complete.map.end$smooth)),size=g_size, colour = "#242B38", cut=TRUE, fill = "#242B38", alpha=tans, adjust = 0.4) + 
  geom_line(aes(RNAmap.incomplete.map.end$end, as.vector(RNAmap.incomplete.map.end$smooth)),size=g_size, colour = "#1D912F", cut=TRUE, fill = "#1D912F", alpha=tans, adjust = 0.4) + 
  ggtitle("eIFA3 read end positions") + 
  xlab("position relative to exon end") + 
  ylab("normalised cDNA counts") + 
  theme(text=element_text(size=6),axis.text=element_text(size=6), axis.title=element_text(size=6,face="plain")) + 
  scale_x_continuous(limits = c(-50, 50))

