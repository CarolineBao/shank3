#this code is used to properly process bed-like annotated GTF files from Gencode to create a gene and feature file to be used in shank3.R
library(dplyr)

FILE_PATH="data/mouse/gencode_mm10/"

genome<-rtracklayer::import('../shank_data/gencode.vM24.annotation.gtf') # read annotated genome GTF file into GRanges

#removing unnecessary columns and repeated rows
genome_df<-as.data.frame(genome, stringsAsFactors= F )
shank_all<-(subset(subset(genome_df, gene_name=="Shank3"), 
                           !(type=="start_codon" | type=="transcript" | type=="CDS" | type=="stop_codon")))
shank_all<-subset(shank_all, select = c(1,2,3,4,5,7,12))
shank_all<-subset(shank_all, !duplicated(shank_all))

#split all_features into gene and features
shank_gene<-subset(shank_all, type=="gene")
shank_features<-data.frame(subset(shank_all, !(type=="gene")), stringsAsFactors = F)

#sorting
shank_features<- shank_features[with(shank_features, order(start)),]

#finding introns 
{
  #add intron as factor
  levels(shank_features$type)<-c(levels(shank_features$type), "intron")

  #if start site doesn't match gene
  if (shank_gene[1,2]<shank_features[1,2]) {
    shank_features <- rbind(shank_features, list(shank_features[1,1], shank_gene[1,2], shank_features[1,2]-1, shank_features[1,2]-shank_gene[1,2], shank_features[1,5], "intron", shank_features[1,7]))
  }
  
  #if end site doesn't match gene
  if (shank_gene[1,3]>shank_features[nrow(shank_features),3]) {
    shank_features <- rbind(shank_features, list(shank_features[1,1], shank_features[nrow(shank_features),3]+1, shank_gene[1,3], shank_gene[1,3]-shank_features[nrow(shank_features),3], shank_features[1,5], "intron", shank_features[1,7]))
  }
  
  #finding introns between features
  for (row in 1:(nrow(shank_features)-1)) {
    #check if there is a gap
    if ((shank_features[row,3] < (shank_features[row+1,2]-1))) {
      shank_features <- rbind(shank_features, list(shank_features[row,1], shank_features[row,3]+1, shank_features[row+1,2]-1, shank_features[row+1,2]-1-shank_features[row,3], shank_features[row,5], "intron", shank_features[row,7]))
    }
  }
  
  #remove duplicates in case UTR and exons lead to repeated introns
  shank_features<- shank_features[with(shank_features, order(start)),]
  shank_features<-subset(shank_features, !duplicated(shank_features))
}

#output
write.table(shank_features, paste(FILE_PATH, "results/shank_features.txt", sep=''), sep = "\t", quote = F, row.names = F)
write.table(shank_gene, paste(FILE_PATH, "results/shank_gene.txt", sep=''), sep = "\t", quote = F, row.names = F)


