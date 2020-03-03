#Uses annotations from UCSC table browser to produce readable files for motif matching. Change gene_features' file path, gene_of_interest, gene_name, and the output file path names to match your needs.
FILE_PATH="data/mouse/gencode_mm9/"

gene_features <- read.table(paste(FILE_PATH,"input_data/shank3_annotated_mm_UCSC.txt", sep=''), header = F, stringsAsFactors = F)
gene_of_interest <-"uc007npi.1"
gene_name<-"Shank3"

#creating the .bed format columns
gene_features<-gene_features %>%  
                      cbind(.,(.$V5-.$V4), rep(gene_name, nrow(gene_features))) %>%
                      subset(str_detect(V10, gene_of_interest)) %>%
                      .[,c(1,4,5,15,7,3,16)]

#rename the column names to match .bed format
colnames(gene_features) <- (c("seqnames", "start", "end", "width", "strand", "type", "gene_name"))

#create dataframe for gene file
gene<-cbind(seqnames=gene_features[1,1], start=gene_features[1,2], end=gene_features[nrow(gene_features),3], 
                width=gene_features[nrow(gene_features),3]-gene_features[1,2]+1, strand=gene_features[1,5], 
                type="gene", gene_name=gene_name)

up_2k_gene <-cbind(seqnames=gene_features[1,1], start=gene_features[1,2]-2000, end=gene_features[nrow(gene_features),3], 
                   width=gene_features[nrow(gene_features),3]-gene_features[1,2]+2001, strand=gene_features[1,5], 
                   type="gene", gene_name=gene_name)

write.table(lapply(as.data.frame(gene_features, stringsAsFactors=F), as.character), paste(FILE_PATH,"input_data/",gene_name,"_features.txt",sep=''), sep = "\t", quote = F, row.names = F)
write.table(lapply(as.data.frame(gene, stringsAsFactors=F), as.character), paste(FILE_PATH,"input_data/",gene_name,"_gene.txt",sep=''), sep = "\t", quote = F, row.names = F)
write.table(lapply(as.data.frame(up_2k_gene, stringsAsFactors=F), as.character), paste(FILE_PATH,"input_data/",gene_name,"_up2kgene.txt",sep=''), sep = "\t", quote = F, row.names = F)

