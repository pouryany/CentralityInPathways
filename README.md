README
================
Pourya Naderi
2/17/2020

# Centrality in Pathways Manual and Walk Through

# Overview

This document contains the walk-through for the accompanying code for
the publication “Revisiting the Use of Graph Centrality Models in
Biological Pathway Analysis”.

# Citation

TBA

# Dependencies and Installation guide

The provided codes have the following dependencies: `graph`, `stringr`,
`dplyr`, `magrittr`,
`CHRONOS`,`pathview`,`KEGGgraph`,`Rgraphviz`,`igraph`,
`org.Mm.eg.db`,`tm`,`stringi`,`stringr`,`dplyr`,`sna`,
`RBGL`,`tidyr`,`AnnotationDbi`,`org.Hs.eg.db`,`annotate`,`biomaRt`,`Hmisc`,`broom`,`xtable`.
Make sure the packages are installed before using. To install the
packages you can use the following code chunk.

``` r
libList <- c("CHRONOS","pathview","KEGGgraph","Rgraphviz","igraph",
             "org.Mm.eg.db","tm","stringi","stringr","dplyr","sna",
             "RBGL","tidyr","AnnotationDbi","org.Hs.eg.db","annotate",
             "biomaRt","Hmisc","broom","xtable")

#source("https://bioconductor.org/biocLite.R")

# for(i in libList) {
#     if(!require(i)){
#         install.packages(i,repos = "http://cran.us.r-project.org")
#     }
#     if(!require(i)){
#         BiocManager::install(i)
#     }
# }
```

## File guides

The folders `human_data` and `mouse_data` contain background data,
preprocessed pathways, and precomputed centrality values. and its
subfolders contains relevant data and code accompanying the package and
the original manuscripts. The folders `HUMAN_PATHWAYLIST` and
`MOUSE_PATHWAYLIST` contain the raw unprocessed XML files of the
original pathways.

The folders `Human Processing` and `Mouse Processing` contain the
scripts. Since the files in the folders are identifical, we will only
review the contents of the `Human Processing` folder.

  - `00-basicFunctions.R`: Contains basic functions as used throughout
    the scripts.
  - `01-00-PathwayParsing.R`: Contains a retrieve script for KEGG
    Pathways and processing them into graph objects. This script
    calculates some basic stats regarding pathways as well.
  - `01-01-CentralityTable.R`: Produces centrality values for all
    qualified pathways using the centrality models, as described in the
    manuscript.  
  - `02-pathwayCleaning.R`: Calculates the quantile score as described
    in the manuscript and filters out the pathways containing less than
    5 cancer genes.
  - `03-1-regressionAnalysis.R`: Computes the regression models for the
    percentage of cancer genes with in each quantile-normalized score
    and generates reports.
  - `03-2-KStests.R`: Kolmogorov-Smirnov test statistics for comparing
    the distribution of cancer genes and non-cancer genes across all
    pathways.  
  - `03-3-pathwayWiseTesting.R`: Contains two-sample test between the
    centrality values of cancer genes and non-cancer genes within each
    pathway.
  - `04-01-pgr_sensitivity.R`: This script generates PageRank centrality
    values across different alpha values.
  - `04-02-pgr_sensitivity_fit.R`: Repeats the regression analysis for
    PageRank centrality values from the last script.

# References

Naderi Yeganeh P, Richardson C, Saule E, Loraine A, Taghi Mostafavi M. Revisiting the use of graph centrality models in biological pathway analysis. BioData Min. 2020;13:5. Published 2020 Jun 16. doi:10.1186/s13040-020-00214-x
