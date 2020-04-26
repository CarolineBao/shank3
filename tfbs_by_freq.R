FILE_PATH <-"data/human/"

shank.match.motif.df<-read.table(paste(FILE_PATH, "results/shank3_region_TF_motifs.txt", sep=''), sep = "\t", header = T, stringsAsFactors = F)

freq <- as.data.frame(table(shank.match.motif.df$motif_nm)) %>%
  .[apply(.!=0, 1, all),] %>%
  cbind(., as.data.frame(shank.match.motif.df$gene_name[1:nrow(.)]))
colnames(freq)<-c("tfbs", "frequency", "gene")

write.table(freq, paste(FILE_PATH, "results/shank3_tfbs_by_freq_hs_ms.txt", sep=''), sep = "\t", quote = F, row.names = F)
