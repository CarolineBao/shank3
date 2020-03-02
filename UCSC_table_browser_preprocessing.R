#uses annotations from UCSC table browser to produce readable files for 
gene_features <- read.table("shank3_annotated_rheMac10_UCSC.txt", header = F, stringsAsFactors = F)
gene_of_interest <-"XM_028827508.1"
gene_name<-"Shank3"
gene_features<-gene_features %>%  
                      cbind(.,(.$V5-.$V4), rep(gene_name, nrow(gene_features))) %>%
                      subset(str_detect(V10, gene_of_interest)) %>%
                      .[,c(1,4,5,15,7,3,16)]

colnames(gene_features) <- (c("seqnames", "start", "end", "width", "strand", "type", "gene_name"))

gene<-cbind(seqnames=gene_features[1,1], start=gene_features[1,2], end=gene_features[nrow(gene_features),3], 
                width=gene_features[nrow(gene_features),3]-gene_features[1,2]+1, strand=gene_features[1,5], 
                type="gene", gene_name=gene_name)

write.table(lapply(as.data.frame(gene_features, stringsAsFactors=F), as.character), paste("output/",gene_name,"_rheMac10_features.txt",sep=''), sep = "\t", quote = F, row.names = F)
write.table(lapply(as.data.frame(gene, stringsAsFactors=F), as.character), paste("output/",gene_name,"_rheMac10_gene.txt",sep=''), sep = "\t", quote = F, row.names = F)

