# Library Imports; Make sure all packages are installed
libList <- c("CHRONOS","pathview","KEGGgraph","Rgraphviz","igraph",
             "org.Mm.eg.db","tm","stringi","stringr","dplyr","sna",
             "RBGL","tidyr","AnnotationDbi","org.Hs.eg.db","annotate",
             "biomaRt","Hmisc","broom","xtable")

source("https://bioconductor.org/biocLite.R")

for(i in libList) {
    if(!require(i)){
        install.packages(i)
    }
    if(!require(i)){
        biocLite(i)
    }
}

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db", version = "3.8")


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

source("00-basicFunctions.R")

# Downloading KEGG PATHWAYS
# Pathwview stores pathways in the working directory

if(!dir.exists("Mouse Processing/MMPATHWAYLIST/")){
    a <- CHRONOS::downloadKEGGPathwayList("mmu")
    pathview::download.kegg(pathway.id = a$Id, species = "mmu",
                            kegg.dir = "./Mouse Processing/MMPATHWAYLIST/")
    paths_list  <- as_data_frame(a)
    write.csv(paths_list, file = "Mouse Processing/pathlist.txt", row.names = F)
}

    a <- read.table("Mouse Processing/pathlist.txt", sep = ",", header = T,
                    colClasses = c("character","character"))
# Retreiving Pathway graphs and information for KGML files
    paths.address  <- paste("Mouse Processing/MMPATHWAYLIST/mmu",a$Id,".xml",sep = "")
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
        lethal.list <- read.csv(file="Mouse Processing/MouseGeneInfo.csv")
        lethal.list <- lethal.list[lethal.list$IMPC_viability_DR7.0 == "Lethal",]
        lethal.list <- lethal.list[,2]
        lethal.list <- as.character(lethal.list)
        mouse.gene.ref <- lethal.list
        mouse.gene.ref <- data_frame(mouse.gene.ref,"Lethal")
        colnames(mouse.gene.ref) <- c("Gene", "Description")




        paths.summary[paths.summary$num_nodes < 21 |paths.summary$num_edges <21,]
        paths.summary[paths.summary$eigen >10,]
        all.path.names     <- pathway.titles
        all.mouse.essential <- data_frame()

j <- 0
k <- 0
for (i in 1:length(graphs.mouse)){


    mouse.Graph  <- graphs.mouse[[i]]
    if(paths.summary  %>% dplyr::select(., num_nodes) %>% slice(.,i) >1000) next
    if(paths.summary  %>% dplyr::select(., num_edges) %>% slice(.,i) >4000) next
    if(paths.summary  %>% dplyr::select(., num_nodes) %>% slice(.,i) <20  ) next
    if(paths.summary  %>% dplyr::select(., num_edges) %>% slice(.,i) <20  ) next
    if(paths.summary  %>% dplyr::select(., eigen)     %>% slice(.,i) > 10 ) next


    node.genes   <- nodes(mouse.Graph)
    b            <- mouse.gene.ref[mouse.gene.ref$Gene %in% node.genes,]
    # List of knockout genes from mouse
    mouse.mat    <- as(mouse.Graph,"matrix")



    #Alpha <- paths.summary  %>% dplyr::select(., eigen)  %>% slice(.,i) %>% pull
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


    #Betweenness Centrality Source
    beet.source.vec     <- sna::betweenness(mouse.mat, cmode = "directed")
    beet.source.rank    <- rank(beet.source.vec, ties.method = "min")
    beet.source.norm    <- zero.one.normalize(beet.source.vec)

    #Betweenness Centrality Sink
    beet.sink.vec     <- sna::betweenness(t(mouse.mat), cmode = "directed")
    beet.sink.rank    <- rank(beet.sink.vec, ties.method = "min")
    beet.sink.norm    <- zero.one.normalize(beet.sink.vec)

    #Betweenness Centrality Source Sink
    beet.ssc.vec     <- beet.source.vec + beet.sink.vec
    beet.ssc.rank    <- rank(beet.ssc.vec, ties.method = "min")
    beet.ssc.norm    <- zero.one.normalize(beet.ssc.vec)


    #Betweenness Centrality Undirected
    beet.und.vec     <- sna::betweenness((mouse.mat), cmode = "undirected")
    beet.und.rank    <- rank(beet.und.vec, ties.method = "min")
    beet.und.norm    <- zero.one.normalize(beet.und.vec)

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
    total.edge   <- pull(paths.summary %>% dplyr::select(., num_edges) %>% slice(.,i))
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

    if(min(beet.source.vec) == max(beet.source.vec)){
        j <- j+1
        next
    }

    mouse.info <- data_frame(pathway.name, total.node,total.edge,node.genes,
                            katz.source.vec, katz.source.rank, katz.source.norm,
                            katz.sink.vec, katz.sink.rank, katz.sink.norm,
                            katz.ssc.vec, katz.ssc.rank, katz.ssc.norm,
                            katz.und.vec, katz.und.rank, katz.und.norm,
                            out.degree, in.degree, all.degree, degree.rank,degree.norm,
                            beet.source.vec, beet.source.rank,beet.source.norm,
                            beet.sink.vec, beet.sink.rank,beet.sink.norm,
                            beet.und.vec, beet.und.rank,beet.und.norm,
                            beet.ssc.vec, beet.ssc.rank, beet.ssc.norm,
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
saveRDS(all.mouse.essential, file = "Mouse Processing/MmsPathwayCentralities.rds")
