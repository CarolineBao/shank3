library(tidyverse)
library(motifmatchr) #bioconductor
library(Matrix)
library(TFBSTools) #bioconductor
library(SummarizedExperiment) #bioconductor
library(BSgenome.Hsapiens.UCSC.hg38) #bioconductor
library(BiocParallel)
library(JASPAR2018)
library(HelloRanges)

register(MulticoreParam(8)) #setting up processing on non-Windows
set.seed(2019) #set random seed
options(warn=-1) #turns off warnings

# method to get JASPAR2018, Getting both human and mouse
#opts <- list()
#opts[["species"]] <- c("Homo sapiens")
jaspar_motifs_hs <- getMatrixSet(JASPAR2018, list("species"="Homo sapiens")) #gets matrix from the named list "species"="Homo sapiens"
jaspar_motifs_ms <- getMatrixSet(JASPAR2018, list("species"="Mus musculus"))
# combining both human and mouse motifs
jaspar_motifs <- c(jaspar_motifs_hs, jaspar_motifs_ms)

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
shank_gene <- read.table("shank3_up2K.txt", header = T, stringsAsFactors = F)
shank_all <- read.table("shank3_allfeatures.txt", header = T, stringsAsFactors = F) #features, aka exons, introns, utr, promoter
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

shank.match.motif.df <- shank.match.motif.ranges[queryHits(feature.overlap)] %>% 
                            as.data.frame %>% 
                            cbind(., as.data.frame(mcols(shank_all[subjectHits(feature.overlap)]))) # add feature info

# join on motif_id to bring in motif meta data, such as species, symbols etc.
shank.match.motif.df <- left_join(shank.match.motif.df, motif_lookup, by = c("group_name" = "motif_id"))

write.table(shank.match.motif.df, "shank3_region_TF_motifs.txt", sep = "\t", quote = F, row.names = F)

#reformatting shank.match.motif.df properly for bedtools
shank.match.motif.df <- makeGRangesFromDataFrame(shank.match.motif.df, keep.extra.columns = T)


#makes windows
modified_makewindows <- function(gene, windowsize) {
    gene_range <- ranges(gene)
    end_value <- (end(gene_range)-start(gene_range))%/%windowsize*windowsize+start(gene_range)+windowsize-1
    gene_temp <- GRanges(seqnames = seqnames(gene), ranges = IRanges(start=start(gene_range), end=end_value, width=end_value-start(gene_range)+1))
    
    output<-tile(gene_temp, width=windowsize)
    output<-as.data.frame(output, stringsAsFactors=F)
    
    output[[5]][[length(output[[1]])]]<-end(gene_range)
    output[[6]][[length(output[[1]])]]<-windowsize-end_value+end(gene_range)
    output<-makeGRangesFromDataFrame(output)
    output
}

intersection_by_bp_window <- function(gene_path, motifs, window_size) {
    
    #Creates windows
    windows <- modified_makewindows(gene, windowsize)
    
    #Calculates frequency of intersection
    ans <- windows
    mcols(ans)$overlap_count <- countOverlaps(windows, motifs, ignore.strand = TRUE)
    write.table(lapply(as.data.frame(ans, stringsAsFactors=F), as.character), 
                paste("shank3_gene_intersect_", window_size, "_bps.txt"), sep="\t", col.names=FALSE, quote = F, row.names = F)
    ans
}

#represents gRanges with count_Overlaps column as a histogram
process_and_graph_overlaps <-function(intersections) {
    to_graph<-as.data.frame(intersections, stringsAsFactors=F)
    to_graph <- subset(to_graph, select = -c(1,2,4,5))
    "plot"(to_graph, type = "histogram")
}

process_and_graph_overlaps(intersection_by_bp_window(shank_gene, shank.match.motif.df, 10))
