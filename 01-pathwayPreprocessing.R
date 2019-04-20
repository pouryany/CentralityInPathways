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
library("org.Hs.eg.db")
library(annotate)
library(biomaRt)
library(Hmisc)
library(broom)
library(xtable)


source("00-basicFunctions.R")

# Downloading KEGG PATHWAYS
# Pathwview stores pathways in the working directory

if(!dir.exists("PATHWAYLIST/")){
    a <- CHRONOS::downloadKEGGPathwayList("hsa")
    pathview::download.kegg(pathway.id = a$Id, kegg.dir = "./PATHWAYLIST/")
    paths_list  <- as_data_frame(a)
    write.csv(paths_list, file = "pathlist.txt", row.names = F)
}

    a <- read.table("pathlist.txt", sep = ",", header = T,
                    colClasses = c("character","character"))
# Retreiving Pathway graphs and information for KGML files
    paths.address  <- paste("PATHWAYLIST/hsa",a$Id,".xml",sep = "")
    graphs.list    <- sapply(paths.address, function(X)(
        KEGGgraph::parseKGML2Graph(X, expandGenes=TRUE)))
    pathway.objs   <- sapply(paths.address, KEGGgraph::parseKGML)
    pathway.titles <- sapply(pathway.objs, getTitle)
    pathway.titles <- unname(pathway.titles)

    names(graphs.list) <- pathway.titles

    graphs.homo    <-  sapply(graphs.list, function(X) ({
        temp <- translateKEGGID2GeneID(nodes(X))
        temp <- getSYMBOL(temp, data='org.Hs.eg')
        nodes(X) <- unname(temp)
        return(X)
    }))



## This line is new and not in submitted version
graphs.homo      <- sapply(graphs.homo, function(X)(RBGL::removeSelfLoops(X)))

### Some basic pathway filterg, leaving out empty graphs
non.empty.homo   <- sapply(graphs.homo, function(X)(length(nodes(X)) != 0))
num.nodes.homo   <- sapply(graphs.homo, function(X)(length(nodes(X))))
non.empty.homo   <- graphs.homo[non.empty.homo]

mtx.collection   <- sapply(non.empty.homo, function(X)(as(X,"matrix")))
num.edges.homo   <- sapply(mtx.collection, sum)
eigen.collection <- lapply(mtx.collection, eigen, only.value = T )

largest.egn.vals <- sapply(eigen.collection,function(X)(largest.eigen(unlist(X))))
egn.vals         <- data_frame(names(largest.egn.vals), largest.egn.vals)
paths.summary    <- data_frame(names(graphs.homo),num.nodes.homo,num.edges.homo)
paths.summary    <- left_join(paths.summary,egn.vals,
                            by = c("names(graphs.homo)" = "names(largest.egn.vals)"))
colnames(paths.summary) <- c("pathway_name", "num_nodes", "num_edges","eigen")


### Generating some pre analysis reports.


        ggplot(paths.summary, aes(x = num_nodes, y = num_edges)) + geom_point()+
            theme_bw() + labs(x= "Nodes", y = "Edges") +
            theme(strip.text = element_text(face="bold", size=20),
                  plot.title = element_text(size = 20),
                  axis.title = element_text(size = 25),
                  legend.text = element_text(size = 20),
                  legend.title=element_blank(),
                  axis.text.y=element_text(size = 20),
                  axis.text.x=element_text(size = 20),
                  axis.ticks.y=element_blank())
        ggsave("images/pathways_node_edge_unfiltered.pdf",
               width = 18, height = 10, units = c("in"))

        ggplot(paths.summary, aes(x = num_nodes, y = eigen)) + geom_point()+
            theme_bw() + labs(x= "Nodes", y = "Eigen-values") +
            theme(strip.text = element_text(face="bold", size=20),
                  plot.title = element_text(size = 20),
                  axis.title = element_text(size = 25),
                  legend.text = element_text(size = 20),
                  legend.title=element_blank(),
                  axis.text.y=element_text(size = 20),
                  axis.text.x=element_text(size = 20),
                  axis.ticks.y=element_blank())
        ggsave("images/pathways_node_eigen_unfiltered.pdf",
               width = 18, height = 10, units = c("in"))


        add.to.row <- list(pos = list(0), command = NULL)
        command <- paste0("\\hline \n \\endhead\n",
                          "\\hline\n",
                          "{\\footnotesize Continued on next page}\n",
                          "\\endfoot\n",
                          "\\endlastfoot\n")
        add.to.row$command <- command
        x.big <- (xtable(paths.summary,digits=c(0,0,0,0,2)))
        print(x.big, hline.after=c(0), add.to.row = add.to.row,
              tabular.environment = "longtable")



# Cancer related genes From Bushman, just an additional option
        cancer.list <- read.table(file = 'data/allOnco_May2018.tsv',
                                  sep = '\t', header = TRUE)
        cancer.list <- as.character(cancer.list$symbol)
        homo.gene.ref <- cancer.list
        homo.gene.ref <- data_frame(homo.gene.ref,"Cancer")
        colnames(homo.gene.ref) <- c("Gene", "Description")


# Cancer related genes From MSigDB
        homo.gene.ref1   <- read.csv("data/msigdblist.csv")
        homo.gene.ref1   <- as_data_frame(homo.gene.ref1)
        homo.gene.ref1$Description <- "Cancer"
        homo.gene.ref1   <-  distinct(homo.gene.ref1,Gene,Description)




#Cancer related genes From Cancer gene consensus
        consensus.list <- read.csv("data/Census_allWed Jun  6 18_56_35 2018.csv")
        consensus.list <- consensus.list[,1]
        homo.gene.ref2 <- data_frame(consensus.list,"Cancer")
        colnames(homo.gene.ref2) <- c("Gene", "Description")


### Selecting union of CGC and MSigDB as the cancer set
        homo.gene.ref <- dplyr::bind_rows(homo.gene.ref1,homo.gene.ref2)
        homo.gene.ref <-  distinct(homo.gene.ref,Gene,Description)

        paths.summary[paths.summary$num_nodes < 21 |paths.summary$num_edges <21,]
        paths.summary[paths.summary$eigen >10,]
        all.path.names     <- pathway.titles
        all.homo.essential <- data_frame()

j <- 0
k <- 0
for (i in 1:length(graphs.homo)){


    homo.Graph  <- graphs.homo[[i]]
    if(paths.summary  %>% dplyr::select(., num_nodes) %>% slice(.,i) >1000) next
    if(paths.summary  %>% dplyr::select(., num_edges) %>% slice(.,i) >4000) next
    if(paths.summary  %>% dplyr::select(., num_nodes) %>% slice(.,i) <20  ) next
    if(paths.summary  %>% dplyr::select(., num_edges) %>% slice(.,i) <20  ) next
    if(paths.summary  %>% dplyr::select(., eigen)     %>% slice(.,i) > 10 ) next


    node.genes   <- nodes(homo.Graph)
    b            <- homo.gene.ref[homo.gene.ref$Gene %in% node.genes,]
    # List of knockout genes from homo
    homo.mat    <- as(homo.Graph,"matrix")



    #Alpha <- paths.summary  %>% dplyr::select(., eigen)  %>% slice(.,i) %>% pull
    #Alpha <- max(min(0.5,(1/Alpha) - 0.1),0.1)
    Alpha <- 0.1

    #Source/Sink Katz Centrality processing, Check the alpha
    cent.mat          <- newpath.centrality(homo.mat,alpha = Alpha, beta = 1)
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
    out.degree   <- rowSums(homo.mat)
    in.degree    <- rowSums(t(homo.mat))
    all.degree   <- in.degree + out.degree
    degree.rank  <- rank(all.degree, ties.method = "min")
    degree.norm  <- zero.one.normalize(all.degree)


    #Betweenness Centrality Source
    beet.source.vec     <- sna::betweenness(homo.mat, cmode = "directed")
    beet.source.rank    <- rank(beet.source.vec, ties.method = "min")
    beet.source.norm    <- zero.one.normalize(beet.source.vec)

    #Betweenness Centrality Sink
    beet.sink.vec     <- sna::betweenness(t(homo.mat), cmode = "directed")
    beet.sink.rank    <- rank(beet.sink.vec, ties.method = "min")
    beet.sink.norm    <- zero.one.normalize(beet.sink.vec)

    #Betweenness Centrality Source Sink
    beet.ssc.vec     <- beet.source.vec + beet.sink.vec
    beet.ssc.rank    <- rank(beet.ssc.vec, ties.method = "min")
    beet.ssc.norm    <- zero.one.normalize(beet.ssc.vec)


    #Betweenness Centrality Undirected
    beet.und.vec     <- sna::betweenness((homo.mat), cmode = "undirected")
    beet.und.rank    <- rank(beet.und.vec, ties.method = "min")
    beet.und.norm    <- zero.one.normalize(beet.und.vec)

    #PageRank Centrality
    igraph.obj   <- igraph::graph_from_adjacency_matrix(homo.mat,mode = "directed")
    igraph.obj2  <- igraph::graph_from_adjacency_matrix(t(homo.mat),mode = "directed")
    igraph.obj3  <- igraph::graph_from_adjacency_matrix(t(homo.mat),mode = "undirected")

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
        source.lap <- semi.laplace(homo.mat)
        sink.lap   <- semi.laplace(t(homo.mat))
        ssc.lap    <- source.lap + sink.lap
        und.lap    <- semi.laplace(igraph::as_adj(igraph.obj3,sparse = F))
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


    # lap.und.vec        <- rowSums(und.lap)
    # lap.und.rank       <- rank(lap.und.vec, ties.method = "min")
    # lap.und.norm       <- zero.one.normalize(lap.und.vec)



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

    # if(as.character(min(lap.und.vec)) == as.character(max(lap.und.vec))){
    #     j <- j+1
    #     next
    # }

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


    homo.info <- data_frame(pathway.name, total.node,total.edge,node.genes,
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


    #homo.info %<>% filter(., node.genes %in% b$Gene) %>%
    #    inner_join(.,b, by = c("node.genes" = "Gene"))

    homo.info <- left_join(homo.info,b, by = c("node.genes" = "Gene"))

    all.homo.essential <- rbind(all.homo.essential, homo.info)
    temp <- all.path.names[[i]]
    print(temp)
}


### Saving and loading made easy
saveRDS(all.homo.essential, file = "pathwayCentralities.rds")
