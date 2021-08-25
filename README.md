The ADARR R package is contains functions for annotating results from the Sprint RNA editing pipeline.

[https://github.com/jumphone/SPRINT](https://github.com/jumphone/SPRINT)

 This annotates each site with the gene identities, genomic regions, and repeat elements overlapping with the sites. This can be integrated with pre-exisiting command line pipelines.  

Using this package requires the GenomicRanges Bioconductor package.  This can be obtained as follows.  

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicRanges")
```
This R package can be installed via Github as follows.

```
install.packages("devtools")
library(devtools)
devtools::install_github("MilesLab/ADARR")
```

Of note, this package requires alignment to the hg38 human genome.  Other annotations may be available in the future. 

Documentation for the functions of ADARR are provided here. 

[https://github.com/MilesLab/ADARR/raw/main/ADARR_0.1.0.pdf](https://github.com/MilesLab/ADARR/raw/main/ADARR_0.1.0.pdf)
