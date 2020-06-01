library(RCy3)

cytoscapePing ()
cytoscapeVersionInfo ()
GENE<-"Shank3"

FILE_PATH<-"/Users/carolinebao/Documents/UROP/Gene\ Therapy/shank3/data/mouse/UCSC_mm10/"
THRESHOLD<-10
freqs <- read.table(paste(FILE_PATH, "results/shank3_tfbs_by_freq_hs_ms_mcm_rr.txt", sep=''), header = T, stringsAsFactors = F)
print(freqs)

#selects specific nodes by their frequency
filtered <- subset(freqs, freqs$frequency > THRESHOLD)  

nodes <- data.frame(id=c(filtered$tfbs, 'Shank3'),
                    stringsAsFactors=FALSE)

edges <- data.frame(source=c(filtered$tfbs),
                    target=c(filtered$gene),
                    weight=c(filtered$frequency), # numeric
                    stringsAsFactors=FALSE)

createNetworkFromDataFrames(nodes, edges, title=paste(GENE, " Frequencies Greater Than", THRESHOLD), collection="DataFrame Example")
