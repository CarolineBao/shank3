# Shank3

Shank3 is a set of R scripts used for identifying transcription factor binding sites along a gene and conducting associated analyses.

## Necessary Installations

```R
install.packages(c("tidyverse", "Matrix"))
chooseCRANmirror() #select mirror number best for you
install.packages("BiocManager")
BiocManager::install(c("motifmatchr", "TFBSTools", "SummarizedExperiment", [GENOME OF INTEREST], "BiocParallel", "JASPAR2020", "HelloRanges")) #replace [GENOME OF INTEREST] with desired genome
```

## Usage

Functions can be run using main.R (identifying transcription factor binding sites from JASPAR2020, finding frequencies for each transcription factor, and binning intersections by window sizes). 

To use main.R:

1) Change header values so that inputs to functions match the information associated with the input data.
2) Go to UCSC Genome Browser to find the position of your desired gene: https://genome.ucsc.edu/
3) Run lines 1-27 of main.R and copy the printed file name.
4) Go to UCSC Table Browser and save the output using the given file name and gene position: https://genome.ucsc.edu/cgi-bin/hgTables
5) Run the rest of main.R to produce motifs and motif intersections.
