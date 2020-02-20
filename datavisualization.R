##???!?
library(tidyverse)


TF.binned.overlap.count <- read.table("intersections.txt", header = F, sep = "\t", dec = ".", stringsAsFactors = F)
TF.binned.overlap.count <- as.data.frame(as.list(TF.binned.overlap.count), stringsAsFactors=F)
TF.binned.overlap.count$V2 <- paste(TF.binned.overlap.count$V2, " to ", TF.binned.overlap.count$V3)
TF.binned.overlap.count <- subset(TF.binned.overlap.count, select = -c(1))
TF.binned.overlap.count <- subset(TF.binned.overlap.count, select = -c(2))
TF.binned.overlap.count$V4<-as.integer(TF.binned.overlap.count$V4)
print(TF.binned.overlap.count)

ggplot(TF.binned.overlap.count, aes(x=TF.binned.overlap.count[2], y=TF.binned.overlap.count[1], fill=TF.binned.overlap.count[2])) + 
  geom_bar(stat="identity", width=.4) +
  geom_text(aes(label=paste0(TF.binned.overlap.count[1], "%\n(", n, ")"), vjust=1.5, colour="white"))

