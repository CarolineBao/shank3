source("scripts/match_tfbs.R")
source("scripts/utils.R")

{
  path <-"/Users/carolinebao/Documents/UROP/Gene\ Therapy/shank3"
  gene <- "vip"
  animal <- "mouse"
  genome <- "BSgenome.Mmusculus.UCSC.mm10"
  motifs <- list("Homo sapiens", "Mus musculus")
  len_upstream=35000
  len_downstream=0
  track="ncbi_refseq"
  feature_source_type="UCSC_table_browser"
  downloads_path="~/Downloads" #Directory files from table browser are downloaded into
}

#initial file must be named as "annotated_[genome_nm]_[track].txt"
setwd(path)
genome_nm<-get_genome_nm(genome)
make_dir_path(gene, animal, genome)

preprocessing_fn<-get_preprocessing_fn(genome, gene, track) #gets the expected file name for the downloaded genome data
print(paste("Expected file name for preprocessing::", preprocessing_fn))
move_download(gene, animal, genome, downloads_path, preprocessing_fn)

match_tfbs(gene, animal, genome, motifs, path, len_upstream, len_downstream, track, feature_source_type)
tfbs_by_freq(gene, animal, genome_nm, len_upstream, len_downstream, path, motifs)

#To-dos: 
# allow for pulling top # of tfbs (eg: )
# allow for pulling specific tfbs data given their names
# add more functionalities (ex: include option to enter own file name, etc)
