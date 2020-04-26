FILE_PATH <-"data/mouse/UCSC_mm10/"

read.table(paste(FILE_PATH, "results/shank3_region_TF_hs_mm_rr_mcm_motifs.txt", sep=''), header = T, stringsAsFactors = F)

freq <- as.data.frame(table(shank.match.motif.df$motif_nm)) %>%
  .[apply(.!=0, 1, all),] %>%
  cbind(., as.data.frame(shank.match.motif.df$gene_name[1:nrow(.)]))
colnames(freq)<-c("tfbs", "frequency", "gene")

write.table(freq, paste(FILE_PATH, "results/shank3_tfbs_by_freq_hs_ms_mcm_rr.txt", sep=''), sep = "\t", quote = F, row.names = F)
