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

Functions can be run using main.R (identifying transcription factor binding sites from JASPAR2020, finding frequencies for each transcription factor, and binning intersections by window sizes of 10, 100, and 500).
