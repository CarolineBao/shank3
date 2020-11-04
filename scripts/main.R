source("scripts/match_tfbs.R")
source("scripts/utils.R")
source("scripts/UCSC_table_browser_preprocessing.R")

{
{
  #Adjust header values to match gene/data information
  path <-"/Users/carolinebao/Documents/Feng UROP (2019-)/Gene Therapy/shank3"
  gene <- "vip_1"
  species <- "mouse"
  genome <- "BSgenome.Mmusculus.UCSC.mm10"
  motifs <- list("Homo sapiens", "Mus musculus")
  len_upstream<-35000
  len_downstream<-0
  track<-"ncbi_refseq"
  feature_source_type<-"UCSC_table_browser"
  downloads_path<-"~/Downloads" #Directory files from table browser are downloaded into
  window_sizes<-list(10, 100, 500)
}

setwd(path)
genome_nm<-get_genome_nm(genome)
make_dir_path(gene, species, genome)

#making/finding initial file and moving it
preprocessing_fn<-get_preprocessing_fn(genome, gene, track) #gets the expected file name for the downloaded genome data to be downloaded from here: https://genome.ucsc.edu/cgi-bin/hgTables
print(paste("Expected file name for preprocessing::", preprocessing_fn))

move_download(gene, species, genome, downloads_path, preprocessing_fn)

#finding tfbs sites and frequencies
match_tfbs(gene, species, genome, motifs, len_upstream, len_downstream, track, feature_source_type)
tfbs_by_freq(gene, species, genome_nm, len_upstream, len_downstream, motifs)

#generating intersections for windows=10, 100, 500
for (window in window_sizes){
  print(paste("Binning frequencies by window size", window))
  intersection_by_bp_window(gene, species, genome, motifs, len_upstream, len_downstream, window_size=window, shift=0, tfbs_name="")
  }
}

#To-dos: 
# allow for pulling top # of tfbs (eg: )
# allow for pulling specific tfbs data given their names
# add more functionalities (ex: include option to enter own file name, etc)

