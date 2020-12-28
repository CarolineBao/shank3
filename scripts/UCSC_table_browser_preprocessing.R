#Uses annotations from UCSC table browser to produce readable files for motif matching. Change gene_features' file path, gene_of_interest, gene_name, and the output file path names to match your needs.
library(tidyverse)
source("scripts/utils.R")

preprocessing <- function (gene, species, genome, len_upstream, len_downstream, track, fn=NULL, write_to_file=TRUE) {
  #Preprocesses data from UCSC table browser

  genome_nm<-get_genome_nm(genome)
  dir_path=paste(paste("data", gene, species, genome_nm, sep='/'), "/", sep="")
  
  if (is.null(fn)){
    #if fn isn't given, use get_preprocessing_fn to get the file name
    fn <- get_preprocessing_fn(genome, gene, track)
  }
  
  #Reads in gene_features and sorts by start
  gene_features <- read_table_helper(fp=paste(dir_path, "input_data/", fn, sep=''))
  gene_features <- gene_features[with(gene_features, order(V5)),]
  
  #Making gene features .bed formatted
  gene_features <- gene_features %>% cbind(.,(.$V5-.$V4), rep(gene, nrow(gene_features))) %>%
    .[,c(1,4,5,15,7,3,16)] %>%
    subset((V3=="exon")) %>%
    subset(!duplicated(gene_features)) %>%
    na.omit
  gene_features <- subset(gene_features, !duplicated(gene_features)) #removing duplicated values in gene_features

  {
    #finding introns between features by checking for gaps between exons
    for (row in 1:(nrow(gene_features)-1)) {
      if (gene_features[row,3] < (gene_features[row+1,2]-1)) {
        #check if there is a gap
        gene_features <- rbind(gene_features, list(gene_features[row,1], gene_features[row,3]+1, gene_features[row+1,2]-1, gene_features[row+1,2]-1-gene_features[row,3], gene_features[row,5], "intron", gene_features[row,7]))
      }
    }
  }
  
  #rename the column names to match .bed format
  colnames(gene_features) <- (c("seqnames", "start", "end", "width", "strand", "type", "gene_name"))
  gene_features<- gene_features[with(gene_features, order(start)),]
  
  #writing out data before flanking information is added
  if (write_to_file){
    fp=paste("data", gene, species, get_genome_nm(genome), "input_data/", sep="/")
    table_writer_checker(paste(fp, gene, "_features_no_flanking.txt", sep=""), gene_features)
    gene_no_flank <-cbind(seqnames=gene_features[1,1], start=gene_features[1,2], end=gene_features[nrow(gene_features),3], 
                          width=gene_features[nrow(gene_features),3]-gene_features[1,2]+1, strand=gene_features[1,5], 
                          type="gene", gene_name=gene)
    table_writer_checker(paste(fp, gene, "_gene_no_flanking.txt", sep=""), gene_no_flank)
  }

  #add in upstream and downstream flanking sections
  if (len_upstream>0){
    gene_features <- rbind(gene_features, list(gene_features[1, 1], gene_features[1,2]-len_upstream, gene_features[1,2]-1, len_upstream, gene_features[1,5], "intron", gene_features[1,7]))
  }

    if (len_downstream>0){
    gene_features <- rbind(gene_features, list(gene_features[1, 1], gene_features[nrow(gene_features),3]+1, gene_features[nrow(gene_features),3]+len_downstream, len_downstream, gene_features[0,5], "3 UTR", gene_features[1,7]))
  }
  gene_features <- gene_features[with(gene_features, order(start)),]
  
  #create dataframe for gene file
  gene<-cbind(seqnames=gene_features[1,1], start=gene_features[1,2], end=gene_features[nrow(gene_features),3], 
              width=gene_features[nrow(gene_features),3]-gene_features[1,2]+1, strand=gene_features[1,5], 
              type="gene", gene_name=gene)
  list(gene=gene, features=gene_features) #returns gene and gene_features as a list
}

