#cran packages
library(tidyverse)
library(Matrix)

#bioconductor packages
library(motifmatchr)
library(TFBSTools)
library(SummarizedExperiment)
library(BiocParallel)
library(JASPAR2020)
library(HelloRanges)

register(MulticoreParam(8)) #setting up processing on non-Windows
set.seed(2019) #set random seed
options(warn=-1) #turns off warnings

make_dir_path <- function(gene, species, genome) {
    #Makes necessary directories for saving and reading data
    {
        cidr <- getwd()
        mkfldr <- paste("data",gene, species, get_genome_nm(genome), "results", "selected_tfbs", "top_tfbs_data", sep="/")
        dir.create(file.path(cidr, mkfldr), recursive = TRUE)
        dir.create(paste("data",gene, species, get_genome_nm(genome), "results", "intersections", sep="/"))
        dir.create(paste("data",gene, species, get_genome_nm(genome), "input_data", sep="/"))
    }
}

match_tfbs <- function (gene, species, genome, motifs, len_upstream, len_downstream, track, feature_source_type) {
    #Inputs
    #   gene (str): gene name
    #   species (str): species name
    #   genome (str): name of genome to be imported
    #   motifs (list): list of species motifs to be imported
    #   path (str): path to the overall folder
    #   len_upstream (int):
    #   len_downstream (int):
    #   feature_source_type (str): where the gene features were sourced from (determines which preprocessing file is run)
    #Output
    #   Saves a file with all of the tfbs locations
    
    #make necessary directories
    make_dir_path(gene, species, genome)
    
    source("scripts/utils.R")
    source(paste("scripts", paste(feature_source_type,"_preprocessing.R", sep=""), sep="/"))
    
    library(genome, character.only = TRUE)
    genome_nm <- get_genome_nm(genome)
    dir_path=paste(paste("data", gene, species, genome_nm, sep='/'), "/", sep="")
    
    #pull all motifs from JASPAR2020
    print("Pulling motifs from JASPAR2020")
    jaspar_motifs <-NULL
    for (a in motifs) {
        if (is.null(jaspar_motifs)) {
            jaspar_motifs<-getMatrixSet(JASPAR2020, list("species"=a))
        }
        else {
            jaspar_motifs <- c(jaspar_motifs, getMatrixSet(JASPAR2020, list("species"=a)))
        }
    }
    
    #represent motifs as string
    motif_nms <- get_motif_nms(motifs)
    
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
    
    # read in gene and all feature annotations from running feature_source_type's respective script
    print("Reading in gene and feature annotations")
    source(paste("scripts", paste(feature_source_type,"_preprocessing.R", sep=""), sep="/"))
    res <- preprocessing(gene, species, genome, len_upstream, len_downstream, track)
    gene_grange <- res$gene
    features_all <- res$features
    gene_grange <- makeGRangesFromDataFrame(as.data.frame(gene_grange), keep.extra.columns = T)
    features_all <- makeGRangesFromDataFrame(as.data.frame(features_all), keep.extra.columns = T)
    
    # match gene with JASPAR2020 motifs
    print("Matching gene with JASPAR2020 motifs and finding overlaps between ranges and annotated features")
    match.motif.pos <- matchMotifs(jaspar_motifs, gene_grange, out = "positions", genome = genome)
    
    # add motif_id to the genomic ranges
    match.motif.ranges <- match.motif.pos %>% as.data.frame %>% 
        dplyr::select(seqnames, start, end, strand, score, group_name) %>%
        makeGRangesFromDataFrame(keep.extra.columns = T)
    
    # get overlaps between motif matching ranges and annotated features
    feature.overlap <- GenomicRanges::findOverlaps(match.motif.ranges, features_all, minoverlap = 10) #use overlaps to find highest frequency
    match.motif.df <- match.motif.ranges[queryHits(feature.overlap)] %>% 
        as.data.frame %>% 
        cbind(., as.data.frame(mcols(features_all[subjectHits(feature.overlap)]))) # add feature info
    
    # join on motif_id to bring in motif meta data, such as species, symbols etc.
    match.motif.df <- left_join(match.motif.df, motif_lookup, by = c("group_name" = "motif_id"))
    print(match.motif.df)
    
    #remove start/end duplicates and saves data
    match.motif.df <-subset(match.motif.df, !duplicated(match.motif.df[,c(2,3,10)]))
    match.motif.df <- apply(match.motif.df,2,as.character) #prevents incompatible characters

    table_writer_checker(paste(dir_path, "results/region_tf_motifs_", motif_nms,"_up_", len_upstream, "_down_", len_downstream,".txt", sep=''), match.motif.df)
    print("Done!")
}

get_motif_nms <- function(motifs) {
    motif_nms <- NULL
    for (motif in motifs) {
        if (is.null(motif_nms)){
            motif_nms <- tolower(substr(motif, 1, 1))
        } else {
            motif_nms <- paste(motif_nms, tolower(substr(motif, 1, 1)), sep="_")
        }
    }
    motif_nms
}

tfbs_by_freq <- function(gene, species, genome_nm, len_upstream, len_downstream, motifs) {
    #Finds the frequencies for each tfbs
    dir_path<-paste("data", gene, species, genome_nm, sep='/')

    motif_nms <- get_motif_nms(motifs)
    fn_descriptor<-paste(motif_nms,"_up_", len_upstream, "_down_", len_downstream, sep="")
    
    match.motif.df<-read_table_helper(paste(dir_path, "/results/region_tf_motifs_",fn_descriptor,".txt", sep=''), sep = "\t", header = T)
    
    freq <- as.data.frame(table(match.motif.df$motif_nm)) %>%
        .[apply(.!=0, 1, all),] %>%
        cbind(., as.data.frame(match.motif.df$gene_name[1:nrow(.)]))
    colnames(freq)<-c("tfbs", "frequency", "gene")
    
    table_writer_checker(paste(dir_path,"/results/tfbs_by_freq_", fn_descriptor,".txt", sep=''), freq)
}


tfbs_by_freq_upstream <- function(gene, species, genome_nm, len_upstream, len_downstream, motifs) {
    #Finds the frequencies for each tfbs
    dir_path<-paste("data", gene, species, genome_nm, sep='/')
    
    motif_nms <- get_motif_nms(motifs)
    fn_descriptor<-paste(motif_nms,"_up_", len_upstream, "_down_", len_downstream, sep="")
    
    get_gene<-read_table_helper(paste(dir_path, "/input_data/",gene,"_gene_no_flanking.txt",sep=''), sep="\t", header=T)
    match.motif.df<-read_table_helper(paste(dir_path, "/results/region_tf_motifs_",fn_descriptor,".txt", sep=''), sep = "\t", header = T)
    match.motif.df<-subset(match.motif.df, match.motif.df$start < get_gene[1, 2])
    
    freq <- as.data.frame(table(match.motif.df$motif_nm)) %>% 
        .[apply(.!=0, 1, all),] %>%
        cbind(., as.data.frame(match.motif.df$gene_name[1:nrow(.)]))
    colnames(freq)<-c("tfbs", "frequency", "gene")
    
    table_writer_checker(paste(dir_path,"/results/tfbs_by_freq_upstream_", fn_descriptor,".txt", sep=''), freq)
}

#makes windows
modified_makewindows <- function(gene, windowsize, shift=0) {
    #Makes bins for a gene and given window size (to be used to count frequencies of tfbs by section on the gene)
    gene_range <- ranges(gene)
    
    #set temporary end value the length of the gene is a multiple of the window size (prevents strange window size creation)
    end_value <- (end(gene_range)-start(gene_range))%/%windowsize*windowsize+start(gene_range)+windowsize-1+shift
    
    #Creates new gene based off of the new end value to be temporarily used to create windows of equal size
    gene_temp <- GRanges(seqnames = seqnames(gene), ranges = IRanges(start=start(gene_range)+shift, end=end_value+shift, width=end_value-start(gene_range)+1))
    
    output<-tile(gene_temp, width=windowsize)    #create windows
    output<-as.data.frame(output, stringsAsFactors=F)    #change windows from dataframe to IRanges
    
    output[1,4] <- start(gene_range) #reset values so they match original ranges
    
    if (end(gene_range)<output[length(output[[1]]),4]){
        #Removes the last window completely if the starting point is after the end of the gene
        output[-nrow(output),]
    }else {
        #Changes the end value on the range of the last window to match the end value of the gene
        output[length(output[[1]]),5]<-end(gene_range)
    }
    
    output[length(output[[1]]),6]<-windowsize-end_value+end(gene_range)
    output<-makeGRangesFromDataFrame(output)
    output
}



intersection_by_bp_window <- function(gene, species, genome, motifs, len_upstream, len_downstream, window_size, shift=0, tfbs_name="") {
    #Bins the intersections by a given window size

    source("scripts/utils.R")
    dir_nm=paste("data", gene, species, get_genome_nm(genome), sep="/")
    descriptor <- paste("_up_", len_upstream, "_down_", len_downstream, sep="")

    #read in and format gene and motif data
    gene_data<-read_table_helper(paste(dir_nm, "input_data", paste(gene, "_gene_no_flanking.txt", sep=''), sep="/"), header=T)
    gene_data[1,2] <- gene_data[1,2]-len_upstream
    gene_data[1,3] <- gene_data[1,3]+len_downstream
    gene_data[1,4] <- gene_data[1,4]+len_upstream+len_downstream
    
    gene_range <- makeGRangesFromDataFrame(gene_data, keep.extra.columns = T)
    motif_df <- as.data.frame(read_table_helper(paste(dir_nm, "/results/region_tf_motifs_",get_motif_nms(motifs), descriptor,".txt", sep=''), sep="\t", header=T))
    
    #Creates windows
    windows <- modified_makewindows(gene_range, window_size, shift)
    
    #Calculates frequency of intersection
    motif_df<- motif_df %>% filter(str_detect(motif_nm, tfbs_name))
    motif_df <- makeGRangesFromDataFrame(motif_df, keep.extra.columns = T)
    ans <- windows
    mcols(ans)$overlap_count <- countOverlaps(windows, motif_df, ignore.strand = TRUE)
    table_writer_checker(paste(dir_nm, "/results/intersections/gene_intersect_frequencies_", get_motif_nms(motifs), descriptor, "_window_", window_size, "_bps",tfbs_name,".txt", sep=""), 
                         lapply(as.data.frame(ans, stringsAsFactors=F), as.character))
    ans
}

intersection_by_bp_window_from_fp <- function(gene_fp, motif_fp, save_fp,len_upstream, len_downstream, window_size, shift=0, tfbs_name="") {
    #Bins the intersections by a given window size
    
    source("scripts/utils.R")
    descriptor <- paste("_up_", len_upstream, "_down_", len_downstream, sep="")
    
    #read in and format gene and motif data
    gene_data<-read_table_helper(paste(dir_nm, "input_data", paste(gene, "_gene_no_flanking.txt", sep=''), sep="/"), header=T)
    gene_data[1,2] <- gene_data[1,2]-len_upstream
    gene_data[1,3] <- gene_data[1,3]+len_downstream
    gene_data[1,4] <- gene_data[1,4]+len_upstream+len_downstream
    
    gene_range <- makeGRangesFromDataFrame(gene_data, keep.extra.columns = T)
    motif_df <- as.data.frame(read_table_helper(motif_fp, sep="\t", header=T))
    
    #Creates windows
    windows <- modified_makewindows(gene_range, window_size, shift)
    
    #Calculates frequency of intersection
    motif_df<- motif_df %>% filter(str_detect(motif_nm, tfbs_name))
    motif_df <- makeGRangesFromDataFrame(motif_df, keep.extra.columns = T)
    ans <- windows
    mcols(ans)$overlap_count <- countOverlaps(windows, motif_df, ignore.strand = TRUE)
    table_writer_checker(save_fp, 
                         lapply(as.data.frame(ans, stringsAsFactors=F), as.character))
    ans
}

process_and_graph_overlaps <-function(intersections) {
    #represents intersections as a histogram
    to_graph<-as.data.frame(intersections, stringsAsFactors=F)
    to_graph <- subset(to_graph, select = -c(1,2,4,5))
    "plot"(to_graph, type = "histogram")
}
