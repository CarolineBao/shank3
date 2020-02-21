
# read GTF file into GRanges
mouse_genome<-rtracklayer::import('gencode.vM24.annotation.gtf')

#removing unnecessary rows and repeated rows
mouse_genome_df<-as.data.frame(mouse_genome, stringsAsFactors= F )
shank_all_mm<-(subset(subset(mouse_genome_df, gene_name=="Shank3"), 
                           !(type=="start_codon" | type=="transcript" | type=="CDS" | type=="stop_codon")))
shank_all_mm<-subset(shank_all_mm, select = c(1,2,3,4,5,7,12))
shank_all_mm<-subset(shank_all_mm, !duplicated(shank_features_mm))

#split all_features into gene and features
shank_gene_mm<-subset(shank_all_mm, type=="gene")
shank_features_mm<-subset(shank_all_mm, !(type=="gene"))

#sorting
shank_features_mm<- shank_features_mm[
  with(shank_features_mm, order(start)),
  ]

#output
write.table(shank_features_mm, "shank_features_mm.txt", sep = "\t", quote = F, row.names = F)
write.table(shank_gene_mm, "shank_gene_mm.txt", sep = "\t", quote = F, row.names = F)
