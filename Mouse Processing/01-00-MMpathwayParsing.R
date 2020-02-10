# Library Imports; Make sure all packages are installed
libList <- c("CHRONOS","pathview","KEGGgraph","Rgraphviz","igraph",
             "org.Mm.eg.db","tm","stringi","stringr","dplyr","sna",
             "RBGL","tidyr","AnnotationDbi","org.Hs.eg.db","annotate",
             "biomaRt","Hmisc","broom","xtable")

# 
# for(i in libList) {
#     if(!require(i)){
#         install.packages(i)
#     }
#     if(!require(i)){
#         BiocManager::install(i)
#     }
# }

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("org.Mm.eg.db", version = "3.8")


library(CHRONOS)
library(pathview)
library(KEGGgraph)
library(Rgraphviz)
library(igraph)
library(org.Mm.eg.db)
library(tm)
library(stringi)
library(stringr)
library(dplyr)
library(sna)
library(RBGL)
library(tidyr)
library("AnnotationDbi")
library("org.Mm.eg.db")
library(annotate)
library(biomaRt)
library(Hmisc)
library(broom)
library(xtable)
library(clusterProfiler)
library(ggplot2)

source("Mouse Processing/00-basicFunctions.R")

# Downloading KEGG PATHWAYS
# Pathwview stores pathways in the working directory

if(!dir.exists("MOUSE_PATHWAYLIST/")){
    a <- CHRONOS::downloadKEGGPathwayList(org = "mmu")
    pathview::download.kegg(pathway.id = a$Id, species = "mmu",
                            kegg.dir = "./MOUSE_PATHWAYLIST/")
    paths_list  <- as_data_frame(a)
    write.csv(paths_list, file = "mouse_data/pathlist.txt", row.names = F)
}

    a <- read.table("mouse_data/pathlist.txt", sep = ",", header = T,
                    colClasses = c("character","character"))
# Retreiving Pathway graphs and information for KGML files
    paths.address  <- paste("MOUSE_PATHWAYLIST/mmu",a$Id,".xml",sep = "")
    graphs.list    <- sapply(paths.address, function(X)(
        KEGGgraph::parseKGML2Graph(X, expandGenes=TRUE)))
    pathway.objs   <- sapply(paths.address, KEGGgraph::parseKGML)
    pathway.titles <- sapply(pathway.objs, getTitle)
    pathway.titles <- unname(pathway.titles)

    names(graphs.list) <- pathway.titles




    graphs.mouse    <-  sapply(graphs.list, function(X) ({
        temp.nodes <- nodes(X)
        temp.nodes <- stringr::str_sub(temp.nodes,start = 5)
        temp <- annotate::getSYMBOL(temp.nodes, data='org.Mm.eg')
        nodes(X) <- unname(temp)
        return(X)
    }))



## This line is new and not in submitted version
graphs.mouse      <- sapply(graphs.mouse, function(X)(RBGL::removeSelfLoops(X)))

### Some basic pathway filterg, leaving out empty graphs
non.empty.mouse   <- sapply(graphs.mouse, function(X)(length(nodes(X)) != 0))
num.nodes.mouse   <- sapply(graphs.mouse, function(X)(length(nodes(X))))
non.empty.mouse   <- graphs.mouse[non.empty.mouse]

mtx.collection   <- sapply(non.empty.mouse, function(X)(as(X,"matrix")))
num.edges.mouse   <- sapply(mtx.collection, sum)
eigen.collection <- lapply(mtx.collection, eigen, only.value = T )

largest.egn.vals <- sapply(eigen.collection,function(X)(largest.eigen(unlist(X))))
egn.vals         <- data_frame(names(largest.egn.vals), largest.egn.vals)
paths.summary    <- data_frame(names(graphs.mouse),num.nodes.mouse,num.edges.mouse)
paths.summary    <- left_join(paths.summary,egn.vals,
                            by = c("names(graphs.mouse)" = "names(largest.egn.vals)"))
colnames(paths.summary) <- c("pathway_name", "num_nodes", "num_edges","eigen")


### Generating some pre analysis reports.

# Mouse Lethal Genes from JAXLABS
        lethal.list <- read.csv(file="mouse_data/MouseGeneInfo.csv")
        lethal.list <- lethal.list[lethal.list$IMPC_viability_DR7.0 == "Lethal",]
        lethal.list <- lethal.list[,2]
        lethal.list <- as.character(lethal.list)
        mouse.gene.ref <- lethal.list
        mouse.gene.ref <- data_frame(mouse.gene.ref,"Lethal")
        colnames(mouse.gene.ref) <- c("Gene", "Description")




        paths.summary[paths.summary$num_nodes < 21 |paths.summary$num_edges <21,]
        paths.summary[paths.summary$eigen >10,]
        all.path.names     <- pathway.titles
        
        
        saveRDS(paths.summary,"mouse_data/paths_summary.RDS")
        saveRDS(all.path.names,"mouse_data/all_path_names.RDS")
        saveRDS(graphs.mouse,"mouse_data/graphs_mouse.RDS")
        saveRDS(mouse.gene.ref,"mouse_data/mouse_gene_ref.RDS")
