library(gplots)
args<-commandArgs(TRUE)
data <- read.table(args[1], header=FALSE, sep=",",stringsAsFactors=F,na.strings=c("NA", "-", "?"))
rnames <- data[,1]
data <- data[,2:250]

data <- as.data.frame(data)

data_matrix <- data.matrix(data)

# following code limits the lowest and highest color to 5%, and 95% of your range, respectively
quantile.range <- quantile(data_matrix, probs = seq(0, 1, 0.01))
palette.breaks <- seq(quantile.range["0%"], quantile.range["99%"], 0.01)

#red
color.palette  <- colorRampPalette(c("#FFFFFF", "#FFBDC0", "#F27F89", "#A12931", "#270004"))(length(palette.breaks) - 1)

pdf(paste(args[2], sep=""))

heatmap.2(
  data_matrix,
  dendrogram = "none",
  scale      = "none",
  trace      = "none",
  Rowv = FALSE,
  Colv = FALSE,
  labRow     = rnames,
  labCol     = c(-50:200),
  col    = color.palette,
  main=paste("",args[3]),
  breaks = palette.breaks,
  key.xlab="pentamer denstiy",
  key.ylab="",
  key.title="",
  cexRow=0.48,
  cexCol=0.48
)
