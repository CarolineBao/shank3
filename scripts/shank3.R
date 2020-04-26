library(tidyverse)
library(motifmatchr) #bioconductor
library(Matrix)
library(TFBSTools) #bioconductor
library(SummarizedExperiment) #bioconductor
library(BSgenome.Hsapiens.UCSC.hg38) #bioconductor
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Mmulatta.UCSC.rheMac10)
library(BiocParallel)
library(JASPAR2018)
library(HelloRanges)

register(MulticoreParam(8)) #setting up processing on non-Windows
set.seed(2019) #set random seed
options(warn=-1) #turns off warnings

FILE_PATH <-"data/human/"

# method to get JASPAR2018, Getting both human and mouse
#opts <- list()
#opts[["species"]] <- c("Homo sapiens")
jaspar_motifs_hs <- getMatrixSet(JASPAR2018, list("species"="Homo sapiens")) #gets matrix from the named list "species"="Homo sapiens"
jaspar_motifs_ms <- getMatrixSet(JASPAR2018, list("species"="Mus musculus"))
jaspar_motifs_mcm <- getMatrixSet(JASPAR2018, list("species"="Macaca mulatta"))
jaspar_motifs_rr <- c(getMatrixSet(JASPAR2018, list("species"="Rattus rattus")), getMatrixSet(JASPAR2018, list("species"="Rattus norvegicus")))
# combining both human and mouse motifs
jaspar_motifs <- c(jaspar_motifs_hs, jaspar_motifs_ms, jaspar_motifs_rr)

# lookup table to join on motif_id to bring together motif information, such as species, symbols, etc.
motif_lookup <- list()
for (m in names(jaspar_motifs)) {
    #motif_lookup[[m]][["ID"]] <- ID(jaspar_motifs[[m]])
    motif_lookup[[m]][["motif_nm"]] <- name(jaspar_motifs[[m]])
    motif_lookup[[m]][["tf_symbol"]] <- ifelse(is.null(tags(jaspar_motifs[[m]])$symbol), "",tags(jaspar_motifs[[m]])$symbol)
    motif_lookup[[m]][["description"]] <- ifelse(is.null(tags(jaspar_motifs[[m]])$description),"", tags(jaspar_motifs[[m]])$description)
    motif_lookup[[m]][["species"]] <- ifelse(is.null(tags(jaspar_motifs[[m]])$species), "", tags(jaspar_motifs[[m]])$species %>% paste0(., collapse = "; "))
}
motif_lookup <- do.call(rbind, motif_lookup)
motif_lookup <- as.data.frame(motif_lookup) %>% rownames_to_column(., "motif_id")

# read in Shank3 gene and all feature annotations
shank_gene <- read.table(paste(FILE_PATH, "input_data/shank_gene.txt", sep=''), header = T, stringsAsFactors = F)
shank_all <- read.table(paste(FILE_PATH, "input_data/shank_features.txt", sep=''), header = T, stringsAsFactors = F) #features, aka exons, introns, utr, promoter
shank_gene <- makeGRangesFromDataFrame(shank_gene, keep.extra.columns = T)
shank_all <- makeGRangesFromDataFrame(shank_all, keep.extra.columns = T)

# match shank gene with jaspar 2018 motifs
shank.match.motif.pos <- matchMotifs(jaspar_motifs, shank_gene, out = "positions", genome = BSgenome.Hsapiens.UCSC.hg38)

# add motif_id to the genomic ranges
shank.match.motif.ranges <- shank.match.motif.pos %>% as.data.frame %>% 
                            dplyr::select(seqnames, start, end, strand, score, group_name) %>%
                            makeGRangesFromDataFrame(keep.extra.columns = T)

# get overlaps between motif matching ranges and shank3 annotated features
feature.overlap <- GenomicRanges::findOverlaps(shank.match.motif.ranges, shank_all, minoverlap = 10) #use overlaps to find highest frequency
print(feature.overlap)
shank.match.motif.df <- shank.match.motif.ranges[queryHits(feature.overlap)] %>% 
                            as.data.frame %>% 
                            cbind(., as.data.frame(mcols(shank_all[subjectHits(feature.overlap)]))) # add feature info

# join on motif_id to bring in motif meta data, such as species, symbols etc.
shank.match.motif.df <- left_join(shank.match.motif.df, motif_lookup, by = c("group_name" = "motif_id"))

#remove start/end duplicates and output data
shank.match.motif.df <-subset(shank.match.motif.df, !duplicated(shank.match.motif.df[,c(2,3,10)]))
write.table(shank.match.motif.df, paste(FILE_PATH, "results/shank3_region_TF_hs_mm_motifs.txt", sep=''), sep = "\t", quote = F, row.names = F)

#makes windows
modified_makewindows <- function(gene, windowsize, shift=0) {
    gene_range <- ranges(gene)
    #set temporary end value so that (length of gene)%windowsize=0
    end_value <- (end(gene_range)-start(gene_range))%/%windowsize*windowsize+start(gene_range)+windowsize-1+shift
    #temporary gene with matching range to temporary end_value
    gene_temp <- GRanges(seqnames = seqnames(gene), ranges = IRanges(start=start(gene_range)+shift, end=end_value+shift, width=end_value-start(gene_range)+1))
    
    output<-tile(gene_temp, width=windowsize)    #create windows
    output<-as.data.frame(output, stringsAsFactors=F)    #change windows from dataframe to IRanges
    
    #reset values so they match original ranges
    output[1,4] <- start(gene_range) 
    if (end(gene_range)<output[length(output[[1]]),4]){
        output[-nrow(output),]
    }else {
        output[length(output[[1]]),5]<-end(gene_range)
    }
    output[length(output[[1]]),6]<-windowsize-end_value+end(gene_range)
    
    output<-makeGRangesFromDataFrame(output)
    output
}

intersection_by_bp_window <- function(gene, motifs, window_size, shift=0, tfbs_name="") {
    #Creates windows
    windows <- modified_makewindows(gene, window_size, shift)
    
    #Calculates frequency of intersection
    motifs<- motifs %>% filter(str_detect(motif_nm, tfbs_name))
    motifs <- makeGRangesFromDataFrame(motifs, keep.extra.columns = T)
    ans <- windows
    mcols(ans)$overlap_count <- countOverlaps(windows, motifs, ignore.strand = TRUE)
    #write.table(lapply(as.data.frame(ans, stringsAsFactors=F), as.character), 
                #paste(FILE_PATH, "results/shank3_gene_intersect_frequencies_", window_size, "_bps",tfbs_name,".txt", sep=''), sep="\t", col.names=FALSE, quote = F, row.names = F)
    print())
    ans
}

#represents gRanges with count_Overlaps column as a histogram
process_and_graph_overlaps <-function(intersections) {
    to_graph<-as.data.frame(intersections, stringsAsFactors=F)
    to_graph <- subset(to_graph, select = -c(1,2,4,5))
    "plot"(to_graph, type = "histogram")
}


results<-intersection_by_bp_window(shank_gene, shank.match.motif.df, 500)

process_and_graph_overlaps(results)
