#Uses annotations from UCSC table browser to produce readable files for motif matching. Change gene_features' file path, gene_of_interest, gene_name, and the output file path names to match your needs.
library(tidyverse)
source("scripts/utils.R")

preprocessing <- function (gene, animal, genome_nm, path, len_upstream, len_downstream, track, fn=NULL, write_to_file=TRUE) {
  #Preprocesses data from UCSC table browser
  #track (str): as shown in UCSC table browser
  genome_nm<-get_genome_nm(genome)
  dir_path=paste(paste(gene, animal, genome_nm, sep='/'), "/", sep="")
  
  if (is.null(fn)){
    fn <- get_preprocessing_fn(genome, gene, track)
  }
  
  gene_features <- tryCatch({read.table(paste(dir_path, "input_data/", fn, sep=''), header = F, stringsAsFactors = F)
  }, error = function(e) {
    stop(paste("Expected ", paste(dir_path, "input_data/", fn, sep=''), " but did not find it. Is the file named correctly?", sep=""))
  })
  
  gene_features <- gene_features[with(gene_features, order(V5)),]
  
  #creating the .bed format columns
  gene_features <- gene_features %>% cbind(.,(.$V5-.$V4), rep(gene, nrow(gene_features))) %>%
    .[,c(1,4,5,15,7,3,16)] %>%
    subset(!(identical(V3, "start_codon") | identical(V3,"transcript") | identical(V3,"CDS") | identical(V3,"stop_codon"))) %>%
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
  
  #writing out data before flanking information is added
  if (write_to_file){
    fp=paste(gene, animal, get_genome_nm(genome), "input_data/", sep="/")
    gene_features<- gene_features[with(gene_features, order(V5)),]
    table_writer_checker(paste(fp, gene, "_features_no_flanking.txt", sep=""), gene_features)
    gene_no_flank <-cbind(seqnames=gene_features[1,1], start=gene_features[1,2], end=gene_features[nrow(gene_features),3], 
                          width=gene_features[nrow(gene_features),3]-gene_features[1,2]+1, strand=gene_features[1,5], 
                          type="gene", gene_name=gene)
    table_writer_checker(paste(fp, gene, "_gene_no_flanking.txt", sep=""), gene_features)
  }
  
  
  #rename the column names to match .bed format
  colnames(gene_features) <- (c("seqnames", "start", "end", "width", "strand", "type", "gene_name"))
  
  #add in upstream and downstream sections
  if (len_upstream>0){
    gene_features <- rbind(gene_features, list(gene_features[0, 1], gene_features[1,2]-len_upstream, gene_features[1,2]-1, len_upstream, gene_features[0,5], "Promoter", gene_features[0,7]))
  }
  
  if (len_downstream>0){
    gene_features <- rbind(gene_features, list(gene_features[0, 1], gene_features[nrow(gene_features),3]+1, gene_features[nrow(gene_features),3]+len_downstream, len_downstream, gene_features[0,5], "3 UTR", gene_features[0,7]))
  }
  
  #create dataframe for gene file
  
  gene<-cbind(seqnames=gene_features[1,1], start=gene_features[1,2], end=gene_features[nrow(gene_features),3], 
              width=gene_features[nrow(gene_features),3]-gene_features[1,2]+1, strand=gene_features[1,5], 
              type="gene", gene_name=gene)
  list(gene=gene, features=gene_features)
}

