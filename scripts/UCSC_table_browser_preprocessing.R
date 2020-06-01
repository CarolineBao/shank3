#Uses annotations from UCSC table browser to produce readable files for motif matching. Change gene_features' file path, gene_of_interest, gene_name, and the output file path names to match your needs.
library(tidyverse)

setwd("/Users/carolinebao/Documents/UROP/Gene Therapy/shank3/")
GENE <- "vip"
FILE_PATH=paste("data_", GENE, "/mouse/", sep='')

gene_features <- read.table(paste(FILE_PATH,"input_data/annotated_UCSC_mm10_gencode_vm_24.txt", sep=''), header = F, stringsAsFactors = F)
gene_of_interest <-"ENSMUST00000045738.4"
gene_features<- gene_features[with(gene_features, order(V5)),]

#creating the .bed format columns
gene_features <- gene_features %>% cbind(.,(.$V5-.$V4), rep(GENE, nrow(gene_features))) %>%
                      #subset(str_detect(V10, gene_of_interest)) %>%
                      .[,c(1,4,5,15,7,3,16)] %>%
                      subset(!(V3=="start_codon" | V3=="transcript" | V3=="CDS" | V3=="stop_codon")) %>%
                      subset(!duplicated(gene_features)) %>%
                      na.omit

gene_features <- subset(gene_features, !duplicated(gene_features))
#finding introns 
{
  #finding introns between features
  for (row in 1:(nrow(gene_features)-1)) {
    #check if there is a gap
    if (gene_features[row,3] < (gene_features[row+1,2]-1)) {
      gene_features <- rbind(gene_features, list(gene_features[row,1], gene_features[row,3]+1, gene_features[row+1,2]-1, gene_features[row+1,2]-1-gene_features[row,3], gene_features[row,5], "intron", gene_features[row,7]))
    }
  }
}
gene_features<- gene_features[with(gene_features, order(V5)),]

#rename the column names to match .bed format
colnames(gene_features) <- (c("seqnames", "start", "end", "width", "strand", "type", "gene_name"))

#create dataframe for gene file
gene<-cbind(seqnames=gene_features[1,1], start=gene_features[1,2], end=gene_features[nrow(gene_features),3], 
                width=gene_features[nrow(gene_features),3]-gene_features[1,2]+1, strand=gene_features[1,5], 
                type="gene", gene_name=GENE)

up_2k_gene <-cbind(seqnames=gene_features[1,1], start=gene_features[1,2]-2000, end=gene_features[nrow(gene_features),3], 
                   width=gene_features[nrow(gene_features),3]-gene_features[1,2]+2001, strand=gene_features[1,5], 
                   type="gene", gene_name=GENE)

write.table(lapply(as.data.frame(gene_features, stringsAsFactors=F), as.character), paste(FILE_PATH,"input_data/",GENE,"_features.txt",sep=''), sep = "\t", quote = F, row.names = F)
write.table(lapply(as.data.frame(gene, stringsAsFactors=F), as.character), paste(FILE_PATH,"input_data/",GENE,"_gene.txt",sep=''), sep = "\t", quote = F, row.names = F)
write.table(lapply(as.data.frame(up_2k_gene, stringsAsFactors=F), as.character), paste(FILE_PATH,"input_data/",GENE,"_up2kgene.txt",sep=''), sep = "\t", quote = F, row.names = F)

