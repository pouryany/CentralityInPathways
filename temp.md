Centrality in Pathways Manual and Walk Through
================
Pourya Naderi Yeganeh
2020-02-17

Overview
========

This document contains the walk-through for the accompanying code for the publication "Revisiting the Use of Graph Centrality Models in Biological Pathway Analysis".

Citation
========

TBA

Dependencies and Installation guide
===================================

The provided codes have the following dependencies:  `graph`, `stringr`, `dplyr`, `magrittr`, `CHRONOS`,`pathview`,`KEGGgraph`,`Rgraphviz`,`igraph`, `org.Mm.eg.db`,`tm`,`stringi`,`stringr`,`dplyr`,`sna`, `RBGL`,`tidyr`,`AnnotationDbi`,`org.Hs.eg.db`,`annotate`,`biomaRt`,`Hmisc`,`broom`,`xtable`. Make sure the packages are installed before using. To install the packages you can use the following code chunk.

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


File guides
-----------

The folders `human_data` and `mouse_data` contain background data, preprocessed pathways, and precomputed centrality values.  and its subfolders contains relevant data and code accompanying the package and the original manuscripts. The folders `HUMAN_PATHWAYLIST` and `MOUSE_PATHWAYLIST` contain  the raw unprocessed XML files of the original pathways.



References
==========
