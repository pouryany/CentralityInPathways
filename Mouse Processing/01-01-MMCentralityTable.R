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


source("Mouse Processing/00-basicFunctions.R")
paths.summary   <- readRDS("mouse_data/paths_summary.RDS")
all.path.names  <- readRDS("mouse_data/all_path_names.RDS")
graphs.mouse    <- readRDS("mouse_data/graphs_mouse.RDS")
mouse.gene.ref  <- readRDS("mouse_data/mouse_gene_ref.RDS")

all.mouse.essential <- data_frame()





all.mouse.essential <- data_frame()




j <- 0
k <- 0
for (i in 1:length(graphs.mouse)){
    
    
    mouse.Graph  <- graphs.mouse[[i]]
    if(paths.summary  %>% dplyr::select(., num_nodes) %>% dplyr::slice(.,i) >1000) next
    if(paths.summary  %>% dplyr::select(., num_edges) %>% dplyr::slice(.,i) >4000) next
    if(paths.summary  %>% dplyr::select(., num_nodes) %>% dplyr::slice(.,i) <20  ) next
    if(paths.summary  %>% dplyr::select(., num_edges) %>% dplyr::slice(.,i) <20  ) next
    if(paths.summary  %>% dplyr::select(., eigen)     %>% dplyr::slice(.,i) > 10 ) next
    
    
    node.genes   <- nodes(mouse.Graph)
    b            <- mouse.gene.ref[mouse.gene.ref$Gene %in% node.genes,]
    # List of knockout genes from mouse
    mouse.mat    <- as(mouse.Graph,"matrix")
    
    
    
    #Alpha <- paths.summary  %>% dplyr::select(., eigen)  %>% dplyr::slice(.,i) %>% pull
    #Alpha <- max(min(0.5,(1/Alpha) - 0.1),0.1)
    Alpha <- 0.1
    
    #Source/Sink Katz Centrality processing, Check the alpha
    cent.mat          <- newpath.centrality(mouse.mat,alpha = Alpha, beta = 1)
    katz.ssc.vec      <- rowSums(cent.mat$SSC)
    katz.ssc.rank     <- rank(katz.ssc.vec, ties.method = "min")
    katz.ssc.norm     <- zero.one.normalize(katz.ssc.vec)
    
    
    #Source Katz Centrality
    katz.source.vec   <- rowSums(cent.mat$Source)
    katz.source.rank  <- rank(katz.source.vec, ties.method = "min")
    katz.source.norm  <- zero.one.normalize(katz.source.vec)
    
    # Sink Katz Centrality
    katz.sink.vec     <- rowSums(cent.mat$Sink)
    katz.sink.rank    <- rank(katz.sink.vec, ties.method = "min")
    katz.sink.norm    <- zero.one.normalize(katz.sink.vec)
    
    #Degree Centrality processing
    out.degree   <- rowSums(mouse.mat)
    in.degree    <- rowSums(t(mouse.mat))
    all.degree   <- in.degree + out.degree
    degree.rank  <- rank(all.degree, ties.method = "min")
    degree.norm  <- zero.one.normalize(all.degree)
    
    
    # #Betweenness Centrality Source
    # cls.source.vec     <- sna::betweenness(mouse.mat, cmode = "directed")
    # cls.source.rank    <- rank(cls.source.vec, ties.method = "min")
    # cls.source.norm    <- zero.one.normalize(cls.source.vec)
    # 
    # #Betweenness Centrality Sink
    # cls.sink.vec     <- sna::betweenness(t(mouse.mat), cmode = "directed")
    # cls.sink.rank    <- rank(cls.sink.vec, ties.method = "min")
    # cls.sink.norm    <- zero.one.normalize(cls.sink.vec)
    # 
    # #Betweenness Centrality Source Sink
    # cls.ssc.vec     <- cls.source.vec + cls.sink.vec
    # cls.ssc.rank    <- rank(cls.ssc.vec, ties.method = "min")
    # cls.ssc.norm    <- zero.one.normalize(cls.ssc.vec)
    # 
    # 
    # #Betweenness Centrality Undirected
    # cls.und.vec     <- sna::betweenness((mouse.mat), cmode = "undirected")
    # cls.und.rank    <- rank(cls.und.vec, ties.method = "min")
    # cls.und.norm    <- zero.one.normalize(cls.und.vec)
    
    #PageRank Centrality
    igraph.obj   <- igraph::graph_from_adjacency_matrix(mouse.mat,mode = "directed")
    igraph.obj2  <- igraph::graph_from_adjacency_matrix(t(mouse.mat),mode = "directed")
    igraph.obj3  <- igraph::graph_from_adjacency_matrix(t(mouse.mat),mode = "undirected")
    
    ### Source Directed pageRank
    pgr.sink.vec   <- igraph::page.rank(igraph.obj,damping = 0.85)
    pgr.sink.vec   <- pgr.sink.vec$vector
    
    ###Source-Sink pageRank
    pgr.source.vec   <- igraph::page.rank(igraph.obj2,damping = 0.85)
    pgr.source.vec   <- pgr.source.vec$vector
    pgr.ssc.vec      <- pgr.source.vec + pgr.sink.vec
    
    ###Undirected pageRank
    pgr.und.vec      <- igraph::page.rank(igraph.obj3,damping = 0.85)
    pgr.und.vec      <- pgr.und.vec$vector
    
    
    pgr.source.rank     <- rank(pgr.source.vec,ties.method = "min")
    pgr.source.norm     <- zero.one.normalize(pgr.source.vec)
    pgr.sink.rank       <- rank(pgr.sink.vec,ties.method = "min")
    pgr.sink.norm       <- zero.one.normalize(pgr.sink.vec)
    pgr.ssc.norm        <- zero.one.normalize(pgr.ssc.vec)
    pgr.ssc.rank        <- rank(pgr.ssc.vec,ties.method = "min")
    pgr.und.norm        <- zero.one.normalize(pgr.und.vec)
    pgr.und.rank        <- rank(pgr.und.vec,ties.method = "min")
    
    #Closeness Centrality Source
    cls.source.vec     <- igraph::closeness(igraph.obj,node.genes, mode = "out")
    cls.source.rank    <- rank(cls.source.vec, ties.method = "min")
    cls.source.norm    <- zero.one.normalize(cls.source.vec)
    
    #Betweenness Centrality Sink
    cls.sink.vec     <- igraph::closeness(igraph.obj,node.genes, mode = "in")
    cls.sink.rank    <- rank(cls.sink.vec, ties.method = "min")
    cls.sink.norm    <- zero.one.normalize(cls.sink.vec)
    
    #Betweenness Centrality Source Sink
    cls.ssc.vec     <- cls.source.vec + cls.sink.vec
    cls.ssc.rank    <- rank(cls.ssc.vec, ties.method = "min")
    cls.ssc.norm    <- zero.one.normalize(cls.ssc.vec)
    
    
    #Betweenness Centrality Undirected
    cls.und.vec     <- igraph::closeness(igraph.obj,node.genes, mode = "all")
    cls.und.rank    <- rank(cls.und.vec, ties.method = "min")
    cls.und.norm    <- zero.one.normalize(cls.und.vec)
    
    
    
    # Undirected katz
    ktz.mat             <-  igraph::as_adj(igraph.obj3,sparse = F)
    katz.und.vec        <-  newpath.centrality(ktz.mat, alpha = 0.1, beta = 0)
    katz.und.vec        <-  rowSums(katz.und.vec$Source)
    katz.und.norm       <-  rank(katz.und.vec, ties.method = "min")
    katz.und.rank       <-  zero.one.normalize(katz.und.vec)
    
    # Semi Laplacian, heat diffusion kernel esque
    tryCatch({
        source.lap <- semi.laplace(mouse.mat)
        sink.lap   <- semi.laplace(t(mouse.mat))
        ssc.lap    <- source.lap + sink.lap
        k <- k+1
    } , error = function(e){print("err"); next}, finally = {print("Grrr")})
    
    
    lap.source.vec     <- rowSums(source.lap)
    lap.source.rank    <- rank(lap.source.vec, ties.method = "min")
    lap.source.norm    <- zero.one.normalize(lap.source.vec)
    
    lap.sink.vec       <- rowSums(sink.lap)
    lap.sink.rank      <- rank(lap.sink.vec, ties.method = "min")
    lap.sink.norm      <- zero.one.normalize(lap.sink.vec)
    
    
    lap.ssc.vec        <- rowSums(ssc.lap)
    lap.ssc.rank       <- rank(lap.ssc.vec, ties.method = "min")
    lap.ssc.norm       <- zero.one.normalize(lap.ssc.vec)
    
    
    
    pathway.name <- all.path.names[[i]]
    pathway.name <- rep(pathway.name,length(node.genes))
    total.node   <- rep(length(node.genes),length(node.genes))
    total.edge   <- pull(paths.summary %>% dplyr::select(., num_edges) %>% dplyr::slice(.,i))
    total.edge   <- rep(total.edge,length(node.genes))
    
    
    # Some quality control
    print(i)
    if(as.character(min(lap.sink.vec)) == as.character(max(lap.sink.vec))){
        print(paste("Catch This Cat",i))
        j <- j+1
        next
    }
    
    if(as.character(min(katz.ssc.vec)) == as.character(max(katz.ssc.vec))){
        j <- j+1
        next
    } else{
        ssc.buckets  <- cut(katz.ssc.vec, breaks = seq(min(katz.ssc.vec),
                                                       max(katz.ssc.vec), len =21),include.lowest = TRUE)
        levels(ssc.buckets) <- 1:20
    }
    
    if(as.character(min(all.degree)) == as.character(max(all.degree))){
        j <- j+1
        next
    } else{deg.buckets  <- cut(all.degree,breaks = seq(min(all.degree), max(all.degree)
                                                       ,len =21), include.lowest = TRUE)
    levels(deg.buckets) <- 1:20
    }
    
    if(min(cls.source.vec) == max(cls.source.vec)){
        j <- j+1
        next
    }
    
    mouse.info <- data_frame(pathway.name, total.node,total.edge,node.genes,
                             katz.source.vec, katz.source.rank, katz.source.norm,
                             katz.sink.vec, katz.sink.rank, katz.sink.norm,
                             katz.ssc.vec, katz.ssc.rank, katz.ssc.norm,
                             katz.und.vec, katz.und.rank, katz.und.norm,
                             out.degree, in.degree, all.degree, degree.rank,degree.norm,
                             cls.source.vec, cls.source.rank,cls.source.norm,
                             cls.sink.vec, cls.sink.rank,cls.sink.norm,
                             cls.und.vec, cls.und.rank,cls.und.norm,
                             cls.ssc.vec, cls.ssc.rank, cls.ssc.norm,
                             pgr.source.vec, pgr.source.rank, pgr.source.norm,
                             pgr.sink.vec, pgr.sink.rank, pgr.sink.norm,
                             pgr.ssc.vec, pgr.ssc.norm, pgr.ssc.rank,
                             pgr.und.vec, pgr.und.norm, pgr.und.rank,
                             lap.source.vec, lap.source.rank, lap.source.norm,
                             lap.sink.vec, lap.sink.rank, lap.sink.norm,
                             lap.ssc.vec, lap.ssc.rank, lap.ssc.norm)
    
    #mouse.info %<>% filter(., node.genes %in% b$Gene) %>%
    #    inner_join(.,b, by = c("node.genes" = "Gene"))
    
    mouse.info <- left_join(mouse.info,b, by = c("node.genes" = "Gene"))
    
    all.mouse.essential <- rbind(all.mouse.essential, mouse.info)
    temp <- all.path.names[[i]]
    print(temp)
}


### Saving and loading made easy
saveRDS(all.mouse.essential, file = "mouse_data/MmsPathwayCentralities.rds")
