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


# Some necessary functions
#
#
#

### Function for calculating Katz-Source-Sink
newpath.centrality  <- function(adj.matrix, alpha, beta){


    eye <- diag(nrow(adj.matrix))
    cent.out  <- solve(eye - alpha * adj.matrix)
    cent.in   <- solve(eye - alpha * t(adj.matrix))
    cent.tot  <- (cent.out) + beta * (cent.in)
    return(list("SSC" = cent.tot, "Sink" =cent.in, "Source" =cent.out))

}

### Given a list of eigenvalues, find the largest real component
largest.eigen       <- function(b) {
    c <-  b %>% unlist()  %>% Im() == 0
    c %<>% subset.default(x = b,.) %>% Re() %>% max()
    return(c)
}

### Normalizes a list of values to zero mean and 1 standard deviation
zero.one.normalize    <- function(cent.list){
    if(max(cent.list) == min(cent.list)){
        return(rep(1,length(cent.list)))
    } else{
    #normalizedVec <- (cent.list - min(cent.list))/(max(cent.list) - min(cent.list))
    normalizedVec <- (cent.list - mean(cent.list))/(sd(cent.list))

    return(normalizedVec)
    }
}

### Semi laplace calculator, works with normalized connectivty matrix
semi.laplace <- function(some.Matrix){

    inv.diag <- 1/(rowSums(some.Matrix) + 0.01)
    inv.diag[!is.finite(inv.diag)] <- 0
    norm.laplace.esque <-  diag(1,length(inv.diag)) - (diag(inv.diag) %*% some.Matrix)
    #norm.laplace.esque <-  Ginv(norm.laplace.esque)
    norm.laplace.esque <-  solve(norm.laplace.esque)
    return(norm.laplace.esque)
}

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
        ggsave("images/pathways_node_eigen_unfiltered.pdf.pdf",
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
cancer.list <- read.table(file = 'data/allOnco_May2018.tsv', sep = '\t', header = TRUE)
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

paths.summary[paths.summary$num_nodes < 21 | paths.summary$num_edges <21,]
paths.summary[paths.summary$eigen >10,]
all.path.names     <- pathway.titles
all.homo.essential <- data_frame()
all.path.names[166]
j <- 0
i <- 10
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



    #Source/Sink Katz Centrality processing
        cent.mat          <- newpath.centrality(homo.mat,alpha = 0.1, beta = 1)
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
        beet.vec     <- sna::betweenness(homo.mat, cmode = "directed")
        beet.rank    <- rank(beet.vec, ties.method = "min")
        beet.norm    <- zero.one.normalize(beet.vec)

    #Betweenness Centrality Sink
        beet.sink.vec     <- sna::betweenness(t(homo.mat), cmode = "directed")
        beet.sink.rank    <- rank(beet.sink.vec, ties.method = "min")
        beet.sink.norm    <- zero.one.normalize(beet.sink.vec)

    #Betweenness Centrality Undirected
        beet.und.vec     <- sna::betweenness((homo.mat), cmode = "undirected")
        beet.und.rank    <- rank(beet.und.vec, ties.method = "min")
        beet.und.norm    <- zero.one.normalize(beet.und.vec)

    #PageRank Centrality
        igraph.obj   <- igraph::graph_from_adjacency_matrix(homo.mat,mode = "directed")
        igraph.obj2  <- igraph::graph_from_adjacency_matrix(t(homo.mat),mode = "directed")
        igraph.obj3  <- igraph::graph_from_adjacency_matrix(t(homo.mat),mode = "undirected")

    ### Source Directed pageRank
        pgr.source.vec   <- igraph::page.rank(igraph.obj,damping = 0.9)
        pgr.source.vec   <- pgr.source.vec$vector

    ###Source-Sink pageRank
        pgr.sink.vec     <- igraph::page.rank(igraph.obj2,damping = 0.9)
        pgr.ssc.vec      <- pgr.source.vec + pgr.sink.vec$vector

    ###Undirected pageRank
        pgr.und.vec      <- igraph::page.rank(igraph.obj3,damping = 0.9)
        pgr.und.vec      <- pgr.und.vec$vector


        pgr.source.rank     <- rank(pgr.source.vec,ties.method = "min")
        pgr.source.norm     <- zero.one.normalize(pgr.source.vec)
        pgr.ssc.norm <- zero.one.normalize(pgr.ssc.vec)
        pgr.dbl.rank <- rank(pgr.ssc.vec,ties.method = "min")
        pgr.und.norm <- zero.one.normalize(pgr.und.vec)
        pgr.und.rank <- rank(pgr.und.vec,ties.method = "min")




    # Semi Laplacian, heat diffusion kernel esque
        tryCatch({
        source.lap <- semi.laplace(homo.mat)
        sink.lap   <- semi.laplace(t(homo.mat))
        ssc.lap    <- source.lap + sink.lap
        und.lap    <- semi.laplace(igraph::as_adj(igraph.obj3,sparse = F))
        k <- k+1
        } , error = function(e){print("err"); next}, finally = {print("Grrr")})


        source.lap.vec     <- rowSums(source.lap)
        source.lap.rank    <- rank(source.lap.vec, ties.method = "min")
        source.lap.norm    <- zero.one.normalize(source.lap.vec)

        sink.lap.vec       <- rowSums(sink.lap)
        sink.lap.rank      <- rank(sink.lap.vec, ties.method = "min")
        sink.lap.norm      <- zero.one.normalize(sink.lap.vec)


        ssc.lap.vec        <- rowSums(ssc.lap)
        ssc.lap.rank       <- rank(ssc.lap.vec, ties.method = "min")
        ssc.lap.norm       <- zero.one.normalize(ssc.lap.vec)


        und.lap.vec        <- rowSums(und.lap)
        ssc.lap.rank       <- rank(und.lap.vec, ties.method = "min")
        ssc.lap.norm       <- zero.one.normalize(und.lap.vec)



    pathway.name <- all.path.names[[i]]
    pathway.name <- rep(pathway.name,length(node.genes))
    total.node   <- rep(length(node.genes),length(node.genes))
    total.edge   <- pull(paths.summary %>% dplyr::select(., num_edges) %>% slice(.,i))
    total.edge   <- rep(total.edge,length(node.genes))


    # Some quality control
    print(i)
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

    if(min(beet.vec) == max(beet.vec)){
        j <- j+1
        next
    } else{
        beet.buckets <- cut(beet.vec,breaks = seq(min(beet.vec), max(beet.vec)
                                                  ,len =21), include.lowest = TRUE)
        levels(beet.buckets) <- 1:20

    }

    homo.path.info <- data_frame(pathway.name, total.node,total.edge,node.genes,
                                 katz.ssc.vec, katz.ssc.rank, katz.ssc.norm,
                                 out.degree, in.degree, all.degree,
                                 degree.rank,degree.norm,
                                 beet.vec, beet.rank,beet.norm,
                                 pgr.source.vec, pgr.source.rank, pgr.source.norm,
                                 pgr.ssc.vec, pgr.ssc.norm, pgr.dbl.rank,
                                 pgr.und.vec, pgr.und.norm, pgr.und.rank,
                                 katz.source.vec, katz.source.rank, katz.source.norm,
                                 ssc.buckets, beet.buckets, deg.buckets)

    #homo.path.info %<>% filter(., node.genes %in% b$Gene) %>%
    #    inner_join(.,b, by = c("node.genes" = "Gene"))

    homo.path.info <- left_join(homo.path.info,b, by = c("node.genes" = "Gene"))

    all.homo.essential <- rbind(all.homo.essential, homo.path.info)
    temp <- all.path.names[[i]]
    print(temp)
}

k
gene.essential <-  all.homo.essential
gene.essential %>% distinct(.,pathway.name)
gene.essential <-  gene.essential %>%  replace_na(list(Description = "Normal"))

# Assigning quantile scores to centrality models
gene.essential %<>% filter(., total.node >20, total.edge > 20) %>%
    mutate(., ssc.quant = round((katz.ssc.rank/total.node)*100),
           deg.quant = round((degree.rank/total.node)*100),
           bet.quant = round((beet.rank/total.node)*100),
           pgr.quant = round((pgr.source.rank/total.node)*100),
           pgr.dbl.quant = round((pgr.dbl.rank/total.node)*100),
           pgr.und.quant = round((pgr.und.rank/total.node)*100),
           ktz.quant = round((katz.source.rank/total.node)*100))

## Total Number of cancer genes found in all genes
##gene.essential %>% filter(., Description == "Cancer") %>% distinct(.,node.genes)



# Pathway-wise checking if the rank of the cancer nodes are different based on centrality
    ### Getting rid of pathways with no cancer-genes
    no.cancers     <- gene.essential %>% group_by(.,pathway.name) %>%
                      filter(.,Description =="Cancer") %>% distinct(., pathway.name)
    gene.essential <- gene.essential %>% filter(., pathway.name %in% unlist(no.cancers))
    no.cancers     <- gene.essential %>% group_by(.,pathway.name) %>%
                      filter(.,Description =="Cancer") %>%
                      summarise(count = n()) %>%
                      filter(count <= 5) %>%
                      distinct(., pathway.name)


    p.vals.process <- gene.essential[!(gene.essential$pathway.name %in% unlist(no.cancers)),]
    distinct(p.vals.process,pathway.name)
    p.vals.process <- p.vals.process %>% mutate(.,pgr.log.vec2 = log(pgr.ssc.vec),
                                                cent.log.vec = log(katz.ssc.vec),
                                                pgr.log.vec = log(pgr.source.vec),
                                                pgr.log.vec3= log(pgr.und.vec))


    # You can change the t.test to wilcox.test in the below formula to compute
    # Nonparametric statistics. In that case, log values would not differ.
    p.vals.process <- p.vals.process %>%
                      gather(., key = "Centrality", value = "cent_value",pgr.log.vec2,
                             pgr.log.vec,pgr.log.vec3,cent.log.vec,pgr.ssc.vec,
                             pgr.source.vec,pgr.und.vec,degree.norm,katz.ssc.vec,katz.source.norm)%>%
                      group_by(pathway.name,Centrality) %>%
                      do(pval =tidy(t.test( cent_value~ Description,
                               alternative = "greater",paired = F,exact=FALSE, data = .))) %>%
                      unnest()  %>%
                      group_by(Centrality) %>%
                      mutate(., fdr = p.adjust(p.value)) %>%
                      filter(., fdr < 0.05)


    # Creating confusion matrix
    mat1 <- table(p.vals.process$Centrality,p.vals.process$pathway.name)
    mat2 <- table(p.vals.process$pathway.name,p.vals.process$Centrality)


    confusion <- as.matrix(mat1) %*% as.matrix(mat2)
    colnames(confusion) <- c("log SS-Katz","SS-Katz","Degree", "Katz","log PageRank","log SS-PageRank",
                             "log Und. PageRank", "PageRank","SS-PageRank", "Und. PageRank")
    rownames(confusion) <- colnames(confusion)

    # The following line only for nonparametrix you have to remove the log values
    # The Wilcox test will not make any difference with the log values
     confusion2 <- confusion[c("SS-Katz","Degree", "Katz", "PageRank",
                              "SS-PageRank", "Und. PageRank"),
                            c("SS-Katz","Degree", "Katz", "PageRank",
                              "SS-PageRank", "Und. PageRank")]

    confusion2[lower.tri(confusion2,diag = F)] <- ""
    confusion[lower.tri(confusion,diag = F)] <- ""
    xtable(confusion2)

    # gene.temp <- gene.essential
    # gene.temp$Description <- factor(gene.temp$Description)
    # levels(gene.temp$Description)
    # gene.temp$Description <- relevel(gene.temp$Description,ref = "Normal")
    # glimpse(gene.essential)
    # aa <- glm(as.factor(Description) ~ (katz.ssc.norm), family = "binomial" ,data = gene.temp)
    # summary(aa)
    # aa$coefficients
    #



    # Analyzing for only pathways with cancer
    gene.essential <- gene.essential[!(gene.essential$pathway.name %in% unlist(no.cancers)),]

    gene.essential %>% distinct(., node.genes)
    gene.essential %>% filter(.,Description == "Cancer") %>% distinct(., node.genes)



    # Plotting Correlations
    biocLite("GGally")
     library(GGally)
    gene.cor.data <- gene.essential %>% dplyr::select(., deg.quant,ssc.quant,
                                                      pgr.quant,pgr.dbl.quant,
                                                      pgr.und.quant)
    ggpairs(gene.cor.data, method = "pearson", columns = 1:ncol(gene.cor.data),
            columnLabels = c("Degree","Source-Sink Katz","PageRank",
                             "Source-Sink PageRank","Undirected PageRank"),
            title = "",
            axisLabels = "internal",
            diag  = list(continuous=wrap("text", labelSize=30)),
            upper = list(continuous=wrap(ggally_cor,size=10)),
            lower = list(continuous = wrap("points", alpha = 0.3,    size=0.1))
            )  + theme_bw()


#####
    #Ignore this block
    {

    abc <-   gene.essential %>% filter(., katz.ssc.vec <20) %>%
        mutate(., katz.ssc.vec = log((katz.ssc.vec)), pgr.ssc.vec = log(pgr.ssc.vec))
    hist(abc$pgr.ssc.vec)
    abc1 <- abc %>%filter(., Description == "Cancer")
    abc2 <- abc %>%filter(., Description == "Normal")

    #abc1 <- mutate(abc1, pgr.source.vec = zero.one.normalize(pgr.source.vec))
    #abc1 <- mutate(abc1, katz.ssc.vec = zero.one.normalize(katz.ssc.vec))

    ks.test(abc1$pgr.ssc.vec,abc2$pgr.ssc.vec, alternative = "less")


    plot(ecdf(abc1$katz.ssc.vec))
    plot(ecdf(abc2$katz.ssc.vec),add = T)

    zz <- cut(abc$pgr.ssc.vec,breaks = 100)
    levels(zz) <- 1:100
    zz2 <- cut(abc$katz.ssc.vec,breaks = 100)
    levels(zz2) <- 1:100

    abc <-  mutate(abc, pgr.quant2 = as.numeric(zz))
    abc <-  mutate(abc, ssc.quant2 = as.numeric(zz2))

    zzz <-  abc %>% filter(., Description =="Normal") %>%
        dplyr::select(.,pgr.quant2) %>% unlist %>% unname %>% ecdf()

    zzz2<- abc %>% filter(., Description =="Cancer") %>%
        dplyr::select(.,pgr.quant2) %>% unlist %>% unname %>% ecdf()

    zzz3 <-  abc %>% filter(., Description =="Normal") %>%
        dplyr::select(.,ssc.quant2) %>% unlist %>% ecdf()
    zzz4 <-  abc %>% filter(., Description =="Cancer") %>%
        dplyr::select(.,ssc.quant2) %>% unlist %>% ecdf()



    plot(zzz)
    plot(zzz2,col = "red", add = T)
    plot(zzz3,col = "blue")
    plot(zzz4,col = "yellow", add = T)

    }
#####


    #### Getting quantiles across the data

    zz <- cut2(gene.essential$pgr.ssc.norm,m = 3, g = 100)
    levels(zz) <- 1:100
    zz2 <- cut2(gene.essential$katz.ssc.norm,m = 3, g = 100)
    levels(zz2) <- 1:100

    zz3 <- cut2(gene.essential$pgr.und.norm,m = 3, g = 100)
    levels(zz3) <- 1:100


    # Creating quantile scores
gene.essential <-  mutate(gene.essential, pgr.quant2 = as.numeric(zz))
gene.essential <-  mutate(gene.essential, ssc.quant2 = as.numeric(zz2))
gene.essential <-  mutate(gene.essential, pgr.quant3 = as.numeric(zz3))



# Comparing the quantile scores
zzz0 <-  gene.essential %>% dplyr::select(.,ssc.quant2) %>%
         unlist %>% unname #%>% ecdf()

zzz1 <-  gene.essential %>% filter(., Description =="Normal") %>%
         dplyr::select(.,ssc.quant2) %>% unlist %>% unname #%>% ecdf()

zzz2<-   gene.essential %>% filter(., Description =="Cancer") %>%
         dplyr::select(.,ssc.quant2) %>% unlist %>% unname #%>% ecdf()

zzz3 <-  gene.essential %>% filter(., Description =="Normal") %>%
         dplyr::select(.,pgr.quant3) %>% unlist  #%>% ecdf()

zzz4 <-  gene.essential %>% filter(., Description =="Cancer") %>%
         dplyr::select(.,pgr.quant3) %>% unlist  #%>% ecdf()


# Just a check between KS Statistics
aa <-ks.test(zzz2,zzz1, alternative = "less")
aa$statistic


# Getting KS Statistics across models
gene.essential2 <- rbind(gene.essential, mutate(gene.essential, Description = "All genes"))

ssk <- ks.test(gene.essential2$ssc.quant2[gene.essential2$Description == "Cancer"],
               gene.essential2$ssc.quant2[gene.essential2$Description == "Normal"],
               alternative = "less")
ssp <- ks.test(gene.essential2$pgr.quant2[gene.essential2$Description == "Cancer"],
               gene.essential2$pgr.quant2[gene.essential2$Description == "Normal"],
               alternative = "less")
udp <- ks.test(gene.essential2$pgr.quant3[gene.essential2$Description == "Cancer"],
               gene.essential2$pgr.quant3[gene.essential2$Description == "Normal"],
               alternative = "less")

ks.pvals <- c(ssp$p.value,udp$p.value,ssk$p.value)
ks.pvals <- formatC(ks.pvals, digits = 3)
ks.pvals <- paste("KS p-value <",ks.pvals)

# Getting plots out!
gene.essential2 <- gene.essential2 %>% gather(.,key = "Type", value = "ranks", pgr.quant2,pgr.quant3,ssc.quant2)
gene.essential2$Type <- as.factor(gene.essential2$Type)
levels(gene.essential2$Type) <- c("Source-Sink PageRank", "Undirected PageRank", "Source-Sink Katz")

ggplot(gene.essential2,aes(ranks,color = Description)) +
    stat_ecdf(geom = "point") + facet_wrap(~Type ,ncol = 3)+theme_bw() + labs(x = "Quantile", y = "Cumulative Density")+
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme(strip.text = element_text(face="bold", size=20),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.title=element_blank(),
          axis.text.y=element_text(size = 20),
          axis.text.x=element_text(size = 20),
          axis.ticks.y=element_blank(),
          legend.position="bottom", legend.box = "horizontal") +
    geom_text(data=data.frame(x=50, y=0, label=ks.pvals,
                Type = c("Source-Sink PageRank", "Undirected PageRank", "Source-Sink Katz")),
                aes(x,y,label=label),size = 7, inherit.aes=FALSE)



ggplot(gene.essential,aes(pgr.und.quant,pgr.quant)) + geom_point()



# Some more plotting

gene.essential3 <- mutate(gene.essential) %>% filter(.,Description == "Cancer") %>% mutate(., Description = "Source-Sink Cancer")
gene.essential4 <- mutate(gene.essential3, pgr.quant2 = pgr.quant3 , Description = "Undirected Cancer")
gene.essential4 <- rbind(gene.essential3, gene.essential4)

ssp.udp <- ks.test(gene.essential4$pgr.quant2[gene.essential4$Description == "Source-Sink Cancer"],
                   gene.essential4$pgr.quant2[gene.essential4$Description == "Undirected Cancer"],
                   alternative = "less")
ssp.udp$p.value

ggplot(gene.essential4,aes(pgr.quant2, color = Description)) +
    stat_ecdf(geom = "point", size =2) + theme_bw() + labs(x = "Quantile", y = "Cumulative Density")+
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme(strip.text = element_text(face="bold", size=20),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.title=element_blank(),
          axis.text.y=element_text(size = 20),
          axis.text.x=element_text(size = 20),
          axis.ticks.y=element_blank(),
          legend.position="bottom", legend.box = "horizontal") +
    annotate("text",x=50, y= 1, size = 10,  label= "KS p-value < 4e-04")

    ks.test(zzz4,zzz2, alternative = "greater")




## The regression analysis

a <- gene.essential %>% group_by(ssc.quant, Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer") %>%
    mutate(., Centrality = "Katz-Source-Sink",quant = ssc.quant)

b <- gene.essential %>% group_by(deg.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer") %>%
    mutate(Centrality = "Degree",quant = deg.quant )

c <- gene.essential %>% group_by(bet.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer") %>%
    mutate(Centrality = "Betweenness" ,quant = bet.quant)

d <- gene.essential %>% group_by(pgr.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer") %>%
    mutate(Centrality = "PageRank",quant = pgr.quant)

e <- gene.essential %>% group_by(ktz.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "Katz",quant = ktz.quant)

f <- gene.essential %>% group_by(pgr.dbl.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "PageRank-Source-Sink",quant = pgr.dbl.quant )

g <- gene.essential %>% group_by(pgr.und.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "PageRank-Undirected",quant = pgr.und.quant )

total <- rbind(a,b,d,e,f,g)



# Reporting

overall.cors <- total %>% group_by(Centrality) %>%
    do(fit = lm(freq ~ as.numeric(quant), data=.)) %>%  glance(fit)

overall.cors <- overall.cors[,c(1,2,3,6)]
overall.cors[,4] <- formatC(unlist(pull(overall.cors[,4])), format = "e", digits = 2)
print(xtable(overall.cors), include.rownames = F)

text.vals <- paste("Adjusted r-squared:",formatC(pull(overall.cors[,2]),digits = 2),
                   "\n p-value:", pull(overall.cors[,4]))

overall.cors$label <- text.vals

ggplot(total, aes(y = freq, x= quant)) + geom_point()+
     geom_smooth(method= "lm") + #geom_smooth(method= "loess", color="green" , fill = "red") +
     facet_wrap(~Centrality ,ncol = 3) +theme_bw()+ labs(x = "Quantile Score", y = "Fraction of Cancer Genes")+
     theme(strip.text = element_text(face="bold", size=20),
          plot.title = element_text(size = 20),
             axis.title = element_text(size = 30),
             legend.text = element_text(size = 9),
             legend.title=element_text(face = "bold", size = 9),
             axis.text.y=element_text(size = 12),
             axis.text.x=element_text(size = 12),
             axis.ticks.y=element_blank()) +  geom_text(data=overall.cors, x = 50, y =0.5,
               aes(x,y,label=label),size = 6, inherit.aes=FALSE)


# Saving namespace

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")


#####
# Ignore this block

?cor.test
total %>% group_by(Centrality) %>% summarise(.,cor.test(freq,as.numeric(quant))$estimate)

h <- abc %>% group_by(pgr.quant2,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>% mutate(Centrality = "Katz Source/Sink")

i <- abc %>% group_by(ssc.quant2,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>% mutate(Centrality = "Katz Source/Sink")

areg <-lm(freq ~ as.numeric(ssc.quant), data = a)
summary(areg)

breg <-lm(freq ~ as.numeric(deg.quant), data = b)
summary(breg)

creg <-lm(freq ~ as.numeric(bet.quant), data = c)
summary(reg)

dreg <-lm(freq ~ as.numeric(pgr.quant), data = d)
summary(dreg)

reg <-lm(freq ~ as.numeric(ktz.quant), data = e)
summary(reg)

reg <-lm(freq ~ as.numeric(pgr.dbl.quant), data = f)
summary(reg)

reg <-lm(freq ~ as.numeric(pgr.und.quant), data = g)
summary(reg)

reg <-lm(freq ~ as.numeric(pgr.quant2), data = h)
summary(reg)

reg <-lm(freq ~ as.numeric(ssc.quant2), data = i)
summary(reg)



plot(a$ssc.quant,a$freq)
plot(b$deg.quant,b$freq)
plot(c$bet.quant,c$freq)
plot(d$pgr.quant,d$freq)
plot(e$ktz.quant,e$freq)
plot(f$pgr.dbl.quant,f$freq)
plot(g$pgr.und.quant,g$freq)
plot(h$pgr.quant2,h$freq)
plot(i$ssc.quant2,i$freq)

sum(b$n )
#####
