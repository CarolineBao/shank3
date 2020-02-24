mouse_genome<-rtracklayer::import('../shank_data/gencode.vM24.annotation.gtf') # read annotated genome GTF file into GRanges

#removing unnecessary rows and repeated rows
mouse_genome_df<-as.data.frame(mouse_genome, stringsAsFactors= F )
shank_all_mm<-(subset(subset(mouse_genome_df, gene_name=="Shank3"), 
                           !(type=="start_codon" | type=="transcript" | type=="CDS" | type=="stop_codon")))
shank_all_mm<-subset(shank_all_mm, select = c(1,2,3,4,5,7,12))
shank_all_mm<-subset(shank_all_mm, !duplicated(shank_all_mm))

#split all_features into gene and features
shank_gene_mm<-subset(shank_all_mm, type=="gene")
shank_features_mm<-data.frame(subset(shank_all_mm, !(type=="gene")), stringsAsFactors = F)

#sorting
shank_features_mm<- shank_features_mm[with(shank_features_mm, order(start)),]

#finding introns 
{
  #add intron as factor
  levels(shank_features_mm$type)<-c(levels(shank_features_mm$type), "intron")

  #if start site doesn't match gene
  if (shank_gene_mm[1,2]<shank_features_mm[1,2]) {
    shank_features_mm <- rbind(shank_features_mm, list(shank_features_mm[1,1], shank_gene_mm[1,2], shank_features_mm[1,2]-1, shank_features_mm[1,2]-shank_gene_mm[1,2], shank_features_mm[1,5], "intron", shank_features_mm[1,7]))
  }
  
  #if end site doesn't match gene
  if (shank_gene_mm[1,3]>shank_features_mm[nrow(shank_features_mm),3]) {
    shank_features_mm <- rbind(shank_features_mm, list(shank_features_mm[1,1], shank_features_mm[nrow(shank_features_mm),3]+1, shank_gene_mm[1,3], shank_gene_mm[1,3]-shank_features_mm[nrow(shank_features_mm),3], shank_features_mm[1,5], "intron", shank_features_mm[1,7]))
  }
  
  #finding introns between features
  for (row in 1:(nrow(shank_features_mm)-1)) {
    #check if there is a gap
    if ((shank_features_mm[row,3] < (shank_features_mm[row+1,2]-1))) {
      shank_features_mm <- rbind(shank_features_mm, list(shank_features_mm[row,1], shank_features_mm[row,3]+1, shank_features_mm[row+1,2]-1, shank_features_mm[row+1,2]-1-shank_features_mm[row,3], shank_features_mm[row,5], "intron", shank_features_mm[row,7]))
    }
  }
  
  #remove duplicates in case UTR and exons lead to repeated introns
  shank_features_mm<- shank_features_mm[with(shank_features_mm, order(start)),]
  shank_features_mm<-subset(shank_features_mm, !duplicated(shank_features_mm))
}

#output
write.table(shank_features_mm, "shank_features_mm.txt", sep = "\t", quote = F, row.names = F)
write.table(shank_gene_mm, "shank_gene_mm.txt", sep = "\t", quote = F, row.names = F)