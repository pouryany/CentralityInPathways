
# We first calculate the PageRank centrality variations for a grid of alpha
# We then do the regression analysis.

rm(list = ls())
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
library(Hmisc)
library(broom)


source("Human Processing/00-basicFunctions.R")
paths.summary  <- readRDS("human_data/paths_summary.RDS")
all.path.names <- readRDS("human_data/all_path_names.RDS")
graphs.hsap    <- readRDS("human_data/graphs_hsap.RDS")
hsap.gene.ref  <- readRDS("human_data/hsap_gene_ref.RDS")

path.items     <- readRDS("human_data/pathItems.RDS")

graphs.hsap    <- graphs.hsap[path.items]
paths.summary  <- paths.summary[paths.summary$pathway_name   %in% path.items,]
all.path.names <- all.path.names[all.path.names %in% path.items]

all.hsap.essential <- data_frame()


j <- 0 
alphas <- seq(0.1,0.95, 0.01)

# Gotta parallelize the nested loops below at some point. 

for (i in 1:length(graphs.hsap)){
    
    hsap.Graph  <- graphs.hsap[[i]]
    if(paths.summary  %>% dplyr::select(., num_nodes) %>% dplyr::slice(.,i) >1000) next
    if(paths.summary  %>% dplyr::select(., num_edges) %>% dplyr::slice(.,i) >4000) next
    if(paths.summary  %>% dplyr::select(., num_nodes) %>% dplyr::slice(.,i) <20  ) next
    if(paths.summary  %>% dplyr::select(., num_edges) %>% dplyr::slice(.,i) <20  ) next
    if(paths.summary  %>% dplyr::select(., eigen)     %>% dplyr::slice(.,i) > 10 ) next
    j <- j+1
    node.genes   <- nodes(hsap.Graph)
    b            <- hsap.gene.ref[hsap.gene.ref$Gene %in% node.genes,]
    # List of knockout genes from hsap
    hsap.mat    <- as(hsap.Graph,"matrix")
    
    
    for(pgr_damp in alphas){
        
        

        #PageRank Centrality
        igraph.obj   <- igraph::graph_from_adjacency_matrix(hsap.mat,mode = "directed")
        igraph.obj2  <- igraph::graph_from_adjacency_matrix(t(hsap.mat),mode = "directed")
        igraph.obj3  <- igraph::graph_from_adjacency_matrix(t(hsap.mat),mode = "undirected")
        
        ### Source Directed pageRank
        pgr.sink.vec   <- igraph::page.rank(igraph.obj,damping = pgr_damp)
        pgr.sink.vec   <- pgr.sink.vec$vector
        
        ###Source-Sink pageRank
        pgr.source.vec   <- igraph::page.rank(igraph.obj2,damping = pgr_damp)
        pgr.source.vec   <- pgr.source.vec$vector
        pgr.ssc.vec      <- pgr.source.vec + pgr.sink.vec
        
        ###Undirected pageRank
        pgr.und.vec      <- igraph::page.rank(igraph.obj3,damping = pgr_damp)
        pgr.und.vec      <- pgr.und.vec$vector
        
        
        pgr.source.rank     <- rank(pgr.source.vec,ties.method = "min")
        pgr.source.norm     <- zero.one.normalize(pgr.source.vec)
        pgr.sink.rank       <- rank(pgr.sink.vec,ties.method = "min")
        pgr.sink.norm       <- zero.one.normalize(pgr.sink.vec)
        pgr.ssc.norm        <- zero.one.normalize(pgr.ssc.vec)
        pgr.ssc.rank        <- rank(pgr.ssc.vec,ties.method = "min")
        pgr.und.norm        <- zero.one.normalize(pgr.und.vec)
        pgr.und.rank        <- rank(pgr.und.vec,ties.method = "min")
        
        
        
        pathway.name <- all.path.names[[i]]
        pathway.name <- rep(pathway.name,length(node.genes))
        total.node   <- rep(length(node.genes),length(node.genes))
        total.edge   <- pull(paths.summary %>% dplyr::select(., num_edges) %>% dplyr::slice(.,i))
        total.edge   <- rep(total.edge,length(node.genes))
        
        
        hsap.info <- data_frame(pathway.name, total.node,total.edge,node.genes,
                                pgr.source.vec, pgr.source.rank, pgr.source.norm,
                                pgr.sink.vec, pgr.sink.rank, pgr.sink.norm,
                                pgr.ssc.vec, pgr.ssc.norm, pgr.ssc.rank,
                                pgr.und.vec, pgr.und.norm, pgr.und.rank, pgr_damp)
        
        
        #hsap.info %<>% filter(., node.genes %in% b$Gene) %>%
        #    inner_join(.,b, by = c("node.genes" = "Gene"))
        
        hsap.info <- left_join(hsap.info,b, by = c("node.genes" = "Gene"))
        
        all.hsap.essential <- rbind(all.hsap.essential, hsap.info)
        temp <- all.path.names[[i]]
        print(temp)
        print(paste0("\n", pgr_damp))
        
    }
    
}

unique(all.hsap.essential$pathway.name)
### Saving and loading made easy
saveRDS(all.hsap.essential, file = "human_data/pgr_sensitivity_analysis.RDS")
