args<-commandArgs(TRUE)

clusters <- read.table(args[1])
clusters$length <- clusters$V3 - clusters$V2
clusters <- clusters[order(clusters$length),]
clusters$V5 <- ""
clusters2 <- clusters[c("V1", "V2", "V3", "length", "V5", "V6")]
write.table(clusters2, args[2], quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
