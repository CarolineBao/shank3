#this code is used to process the downloaded table from Ensembl to create a gene and feature file to be used in shank3.R

library(dplyr)
library(stringr)
shank_all_features<-read.table("original_macaca_mulatta.csv", header = T, stringsAsFactors = F, sep=",") # read csv file

#removing unnecessary columns
shank_all_features <- shank_all_features[, c(3,4,7,2)]

for (row in 1:(nrow(shank_all_features))) {
  #don't know why it's skipping over first exon
  if (is.na(shank_all_features[row,1]) | (shank_all_features[row,1]=="")) {
    #removing 5' and 3' UTR's
    shank_all_features <- shank_all_features[-row,]
  }
  else{
    
    if (str_detect(shank_all_features[row, 4], "Intron")) {
      shank_all_features[row,4]<-"intron"
    }
    else{
      shank_all_features[row,4]<-"exon"
    }
  }
}
#manually replace exon name
shank_all_features[1,4] <-"exon"

# add in chr, strand, type and rename headers(located on chr 10)
if (shank_all_features[1,1]>shank_all_features[1,2]) {
  strand <- '-'
}else{
  strand<-'+'
}
seqnames <- rep('chr10', nrow(shank_all_features))
strand <- rep(strand, nrow(shank_all_features))
gene_name<-rep('Shank3', nrow(shank_all_features))
shank_all_features <- cbind(shank_all_features, seq_names, strand, gene_name)[, c(5,2,1,3,6,4,7)] #located on negative strand though (flipped start+end values)
colnames(shank_all_features) <- (c("seqnames", "start", "end", "width", "strand", "type", "gene_name"))

#reformat start and end
shank_all_features[,2] <- as.numeric(gsub(",", "", shank_all_features[,2]))
shank_all_features[,3] <- as.numeric(gsub(",", "", shank_all_features[,3]))

#create gene dataframe manually (44 TO 48 DOESN'T WORK RIGHT NOW)
if (strand[1]=='-'){
  shank_gene <- cbind("seqnames"=shank_all_features[1,1], "start"=shank_all_features[nrow(shank_all_features), 3], "end"=shank_all_features[1, 2], 
                     "width"=(shank_all_features[1, 2]-shank_all_features[nrow(shank_all_features), 3]+1), "strand"='-', 
                      "type"="gene", "gene_name"="Shank3")
}else{
  shank_gene <- cbind("seqnames"=shank_all_features[1,1], "start"=shank_all_features[nrow(shank_all_features), 3], "end"=shank_all_features[1, 2],  
                      "width"=(shank_all_features[nrow(shank_all_features), 3]-shank_all_features[1, 2]+1), "strand"='+', 
                      "type"="gene", "gene_name"="Shank3")
}

#output
write.table(shank_all_features, "shank_ensembl_features.txt", sep = "\t", quote = F, row.names = F)
write.table(shank_gene, "shank_ensembl_gene.txt", sep = "\t", quote = F, row.names = F)
