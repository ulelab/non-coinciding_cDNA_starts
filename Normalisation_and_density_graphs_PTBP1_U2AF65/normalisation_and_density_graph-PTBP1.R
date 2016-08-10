library("ggplot2")
library("smoother")

args<-commandArgs(TRUE)
smoothing_window <- 10

##################################################
# read start positions at the start of the cluster

# importing the data
PTB.start.1 <- read.table(args[1], sep='\t')
PTB.start.1$V4 <- 1   #each read count as 1

# getting the maximum peak for each transcript
PTB.start.1.sum <- aggregate(PTB.start.1$V4, list(PTB.start.1$V1, PTB.start.1$V2, PTB.start.1$V6, PTB.start.1$V12), FUN=sum, na.rm=TRUE) # for each reference we get crosslinks summed together
PTB.start.1.max <- aggregate(PTB.start.1.sum$x, list(PTB.start.1.sum$Group.4,PTB.start.1.sum$Group.3), FUN=max, na.rm=TRUE) #for each reference we get maximum

#normalization by the maximum peak and sum of all cDNAs 
colnames(PTB.start.1.max)[1] <- "Group.4"
PTB.start.1.norm <- merge(PTB.start.1.sum, PTB.start.1.max, by="Group.4") #we merge sum with maximum
PTB.start.1.norm$norm <- (PTB.start.1.norm$x.x / PTB.start.1.norm$x.y) / sum(PTB.start.1$V4) #normalize each sum with maximum
PTB.start.1.norm$map <- ifelse(PTB.start.1.norm$Group.3 == '-', PTB.start.1.norm$Group.4 - PTB.start.1.norm$Group.2.x , PTB.start.1.norm$Group.2.x - PTB.start.1.norm$Group.4)  #for minus strand we need to switch positions '+':xnt - BP, '-': BP-xnt
PTB.start.1.map <- aggregate(norm ~ map, data=PTB.start.1.norm, FUN=sum) #we sum together cDNA counts for each position

# normalisation by the average cDNAs in region 25 nt downstream from cluster end
PTB.start.1.map.intronic <- PTB.start.1.map[which(PTB.start.1.map$map > 25),]
PTB.start.1.map$norm <- PTB.start.1.map$norm / mean(PTB.start.1.map.intronic$norm)  #normalize normlaized values with and average count downstream from Low Complexity region

# smoothing
PTB.start.1.map$smooth <- smth(PTB.start.1.map$norm, window = smoothing_window, method = "gaussian")


##################################################
# read end positions at the start of the cluster

# importing the data
PTB.start.2 <- read.table(args[1], sep='\t')
PTB.start.2$V4 <- 1   #each read count as 1

# getting the maximum peak for each transcript
PTB.start.2.sum <- aggregate(PTB.start.2$V4, list(PTB.start.2$V1, PTB.start.2$V2, PTB.start.2$V6, PTB.start.2$V12), FUN=sum, na.rm=TRUE) # for each reference we get crosslinks summed together
PTB.start.2.max <- aggregate(PTB.start.2.sum$x, list(PTB.start.2.sum$Group.4,PTB.start.2.sum$Group.3), FUN=max, na.rm=TRUE) #for each reference we get maximum

#normalization by the maximum peak and sum of all cDNAs 
colnames(PTB.start.2.max)[1] <- "Group.4"
PTB.start.2.norm <- merge(PTB.start.2.sum, PTB.start.2.max, by="Group.4") #we merge sum with maximum
PTB.start.2.norm$norm <- (PTB.start.2.norm$x.x / PTB.start.2.norm$x.y) / sum(PTB.start.2$V4) #normalize each sum with maximum
PTB.start.2.norm$map <- ifelse(PTB.start.2.norm$Group.3 == '-', PTB.start.2.norm$Group.4 - PTB.start.2.norm$Group.2.x , PTB.start.2.norm$Group.2.x - PTB.start.2.norm$Group.4)  #for minus strand we need to switch positions '+':xnt - BP, '-': BP-xnt
PTB.start.2.map <- aggregate(norm ~ map, data=PTB.start.2.norm, FUN=sum) #we sum together cDNA counts for each position

# normalisation by the average cDNAs in region 25 nt downstream from cluster end
PTB.start.2.map.intronic <- PTB.start.2.map[which(PTB.start.2.map$map > 25),]
PTB.start.2.map$norm <- PTB.start.2.map$norm / mean(PTB.start.2.map.intronic$norm)  #normalize normlaized values with and average count downstream from Low Complexity region

# smoothing
PTB.start.2.map$smooth <- smth(PTB.start.2.map$norm, window = smoothing_window, method = "gaussian")


##################################################
# read start positions at the end of the cluster

# importing the data
PTB.start.3 <- read.table(args[1], sep='\t')
PTB.start.3$V4 <- 1   #each read count as 1

# getting the maximum peak for each transcript
PTB.start.3.sum <- aggregate(PTB.start.3$V4, list(PTB.start.3$V1, PTB.start.3$V2, PTB.start.3$V6, PTB.start.3$V12), FUN=sum, na.rm=TRUE) # for each reference we get crosslinks summed together
PTB.start.3.max <- aggregate(PTB.start.3.sum$x, list(PTB.start.3.sum$Group.4,PTB.start.3.sum$Group.3), FUN=max, na.rm=TRUE) #for each reference we get maximum

#normalization by the maximum peak and sum of all cDNAs 
colnames(PTB.start.3.max)[1] <- "Group.4"
PTB.start.3.norm <- merge(PTB.start.3.sum, PTB.start.3.max, by="Group.4") #we merge sum with maximum
PTB.start.3.norm$norm <- (PTB.start.3.norm$x.x / PTB.start.3.norm$x.y) / sum(PTB.start.3$V4) #normalize each sum with maximum
PTB.start.3.norm$map <- ifelse(PTB.start.3.norm$Group.3 == '-', PTB.start.3.norm$Group.4 - PTB.start.3.norm$Group.2.x , PTB.start.3.norm$Group.2.x - PTB.start.3.norm$Group.4)  #for minus strand we need to switch positions '+':xnt - BP, '-': BP-xnt
PTB.start.3.map <- aggregate(norm ~ map, data=PTB.start.3.norm, FUN=sum) #we sum together cDNA counts for each position

# normalisation by the average cDNAs in region 25 nt downstream from cluster end
PTB.start.3.map.intronic <- PTB.start.3.map[which(PTB.start.3.map$map > 25),]
PTB.start.3.map$norm <- PTB.start.3.map$norm / mean(PTB.start.3.map.intronic$norm)  #normalize normlaized values with and average count downstream from Low Complexity region

# smoothing
PTB.start.3.map$smooth <- smth(PTB.start.3.map$norm, window = smoothing_window, method = "gaussian")


##################################################
# read end positions at the end of the cluster

# importing the data
PTB.start.4 <- read.table(args[1], sep='\t')
PTB.start.4$V4 <- 1   #each read count as 1

# getting the maximum peak for each transcript
PTB.start.4.sum <- aggregate(PTB.start.4$V4, list(PTB.start.4$V1, PTB.start.4$V2, PTB.start.4$V6, PTB.start.4$V12), FUN=sum, na.rm=TRUE) # for each reference we get crosslinks summed together
PTB.start.4.max <- aggregate(PTB.start.4.sum$x, list(PTB.start.4.sum$Group.4,PTB.start.4.sum$Group.3), FUN=max, na.rm=TRUE) #for each reference we get maximum

#normalization by the maximum peak and sum of all cDNAs 
colnames(PTB.start.4.max)[1] <- "Group.4"
PTB.start.4.norm <- merge(PTB.start.4.sum, PTB.start.4.max, by="Group.4") #we merge sum with maximum
PTB.start.4.norm$norm <- (PTB.start.4.norm$x.x / PTB.start.4.norm$x.y) / sum(PTB.start.4$V4) #normalize each sum with maximum
PTB.start.4.norm$map <- ifelse(PTB.start.4.norm$Group.3 == '-', PTB.start.4.norm$Group.4 - PTB.start.4.norm$Group.2.x , PTB.start.4.norm$Group.2.x - PTB.start.4.norm$Group.4)  #for minus strand we need to switch positions '+':xnt - BP, '-': BP-xnt
PTB.start.4.map <- aggregate(norm ~ map, data=PTB.start.4.norm, FUN=sum) #we sum together cDNA counts for each position

# normalisation by the average cDNAs in region 25 nt downstream from cluster end
PTB.start.4.map.intronic <- PTB.start.4.map[which(PTB.start.4.map$map > 25),]
PTB.start.4.map$norm <- PTB.start.4.map$norm / mean(PTB.start.4.map.intronic$norm)  #normalize normlaized values with and average count downstream from Low Complexity region

# smoothing
PTB.start.4.map$smooth <- smth(PTB.start.4.map$norm, window = smoothing_window, method = "gaussian")


# draw all figures 
tans <- 1.0
g_size <- 0.4
adj <- 0.4
cluster.start1.gg <- ggplot() + theme_bw() + 
  geom_line(aes(PTB.start.1.map$map, as.vector(PTB.start.1.map$smooth)),size=g_size, colour =  "#242B38", cut=TRUE, fill =  "#242B38", alpha=tans, adjust = 0.4) + 
  ggtitle("PTB-iCLIP read starts") + 
  xlab("position relative cluster start") + 
  ylab("normalized density of cDNA counts") + 
  theme(text=element_text(size=6),axis.text=element_text(size=6), axis.title=element_text(size=6,face="plain")) + 
  scale_x_continuous(limits = c(-75, 75)) 
cluster.start1.gg

cluster.start2.gg <- ggplot() + theme_bw() + 
  geom_line(aes(PTB.start.2.map$map, as.vector(PTB.start.2.map$smooth)),size=g_size, colour =  "#242B38", cut=TRUE, fill =  "#242B38", alpha=tans, adjust = 0.4) + 
  ggtitle("PTB-iCLIP read ends") + 
  xlab("position relative cluster start") + 
  ylab("normalized density of cDNA counts") + 
  theme(text=element_text(size=6),axis.text=element_text(size=6), axis.title=element_text(size=6,face="plain")) + 
  scale_x_continuous(limits = c(-75, 75))
cluster.start2.gg

cluster.end1.gg <- ggplot() + theme_bw() + 
  geom_line(aes(PTB.start.3.map$map, as.vector(PTB.start.3.map$smooth)),size=g_size, colour =  "#242B38", cut=TRUE, fill =  "#242B38", alpha=tans, adjust = 0.4) + 
  ggtitle("PTB-iCLIP read starts") + 
  xlab("position relative cluster end") + 
  ylab("normalized density of cDNA counts") + 
  theme(text=element_text(size=6),axis.text=element_text(size=6), axis.title=element_text(size=6,face="plain")) + 
  scale_x_continuous(limits = c(-75, 75))
cluster.end1.gg

cluster.end2.gg <- ggplot() + theme_bw() + 
  geom_line(aes(PTB.start.4.map$map, as.vector(PTB.start.4.map$smooth)),size=g_size, colour =  "#242B38", cut=TRUE, fill =  "#242B38", alpha=tans, adjust = 0.4) + 
  ggtitle("PTB-iCLIP read ends") + 
  xlab("position relative cluster end") + 
  ylab("normalized density of cDNA counts") + 
  theme(text=element_text(size=6),axis.text=element_text(size=6), axis.title=element_text(size=6,face="plain")) + 
  scale_x_continuous(limits = c(-75, 75)) 
cluster.end2.gg


