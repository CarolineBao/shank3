source("scripts/match_tfbs.R")

path <-"/Users/carolinebao/Documents/tfbs_finder"
gene <- "vip"
animal <- "mouse"
genome <- "BSgenome.Mmusculus.UCSC.mm10"
motifs <- list("Homo sapiens", "Mus musculus")
len_upstream=0
len_downstream=0
track="gencode_vm_24"
feature_source_type="UCSC_table_browser"

#initial file must be named as "annotated_[genome_nm]_[track].txt"

setwd(path)

genome_nm<-get_genome_nm(genome)

make_dir_path(gene, animal, genome)
match_tfbs(gene, animal, genome, motifs, path, len_upstream, len_downstream, track, feature_source_type)
tfbs_by_freq(gene, animal, genome_nm, path, motifs)


