library(tidyverse)
setwd("/Users/carolinebao/Documents/UROP/Gene Therapy/shank3/")
GENE<-"pvalb"
FILE_PATH <-paste("data_",GENE,"/mouse/", sep='')


match.motif.df<-read.table(paste(FILE_PATH, "results/regions_TF_motifs_hs_ms.txt", sep=''), sep = "\t", header = T, stringsAsFactors = F)

freq <- as.data.frame(table(match.motif.df$motif_nm)) %>%
  .[apply(.!=0, 1, all),] %>%
  cbind(., as.data.frame(match.motif.df$gene_name[1:nrow(.)]))
colnames(freq)<-c("tfbs", "frequency", "gene")

write.table(freq, paste(FILE_PATH, "results/tfbs_by_freq_hs_ms.txt", sep=''), sep = "\t", quote = F, row.names = F)
