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


source("Human Processing/00-basicFunctions.R")

# Downloading KEGG PATHWAYS
# Pathwview stores pathways in the working directory

if(!dir.exists("PATHWAYLIST/")){
    a <- CHRONOS::downloadKEGGPathwayList("hsa")
    pathview::download.kegg(pathway.id = a$Id, kegg.dir = "./PATHWAYLIST/")
    paths_list  <- as_data_frame(a)
    write.csv(paths_list, file = "human_data/pathlist.txt", row.names = F)
}

a <- read.table("human_data/pathlist.txt", sep = ",", header = T,
                colClasses = c("character","character"))
# Retreiving Pathway graphs and information for KGML files
paths.address  <- paste("PATHWAYLIST/hsa",a$Id,".xml",sep = "")
graphs.list    <- sapply(paths.address, function(X)(
    KEGGgraph::parseKGML2Graph(X, expandGenes=TRUE)))
pathway.objs   <- sapply(paths.address, KEGGgraph::parseKGML)
pathway.titles <- sapply(pathway.objs, getTitle)
pathway.titles <- unname(pathway.titles)

names(graphs.list) <- pathway.titles

graphs.hsap    <-  sapply(graphs.list, function(X) ({
    temp <- translateKEGGID2GeneID(nodes(X))
    temp <- getSYMBOL(temp, data='org.Hs.eg')
    nodes(X) <- unname(temp)
    return(X)
}))



## This line is new and not in submitted version
graphs.hsap      <- sapply(graphs.hsap, function(X)(RBGL::removeSelfLoops(X)))

### Some basic pathway filterg, leaving out empty graphs
length(graphs.hsap)
non.empty.hsap   <- sapply(graphs.hsap, function(X)(length(nodes(X)) != 0))
num.nodes.hsap   <- sapply(graphs.hsap, function(X)(length(nodes(X))))
non.empty.hsap   <- graphs.hsap[non.empty.hsap]
length(non.empty.hsap)

mtx.collection   <- sapply(non.empty.hsap, function(X)(as(X,"matrix")))
num.edges.hsap   <- sapply(mtx.collection, sum)
eigen.collection <- lapply(mtx.collection, eigen, only.value = T )

largest.egn.vals <- sapply(eigen.collection,function(X)(largest.eigen(unlist(X))))
egn.vals         <- data_frame(names(largest.egn.vals), largest.egn.vals)
paths.summary    <- data_frame(names(graphs.hsap),num.nodes.hsap,num.edges.hsap)
paths.summary    <- left_join(paths.summary,egn.vals,
                              by = c("names(graphs.hsap)" = "names(largest.egn.vals)"))
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
hsap.gene.ref <- cancer.list
hsap.gene.ref <- data_frame(hsap.gene.ref,"Cancer")
colnames(hsap.gene.ref) <- c("Gene", "Description")


# Cancer related genes From MSigDB
hsap.gene.ref1   <- read.csv("human_data/msigdblist.csv")
hsap.gene.ref1   <- as_tibble(hsap.gene.ref1)
hsap.gene.ref1$Description <- "Cancer"
hsap.gene.ref1   <-  distinct(hsap.gene.ref1,Gene,Description)




#Cancer related genes From Cancer gene consensus
consensus.list <- read.csv("human_data/Census_allWed Jun  6 18_56_35 2018.csv")
consensus.list <- consensus.list[,1]
hsap.gene.ref2 <- data_frame(consensus.list,"Cancer")
colnames(hsap.gene.ref2) <- c("Gene", "Description")


### Selecting union of CGC and MSigDB as the cancer set
hsap.gene.ref <- dplyr::bind_rows(hsap.gene.ref1,hsap.gene.ref2)
hsap.gene.ref <-  distinct(hsap.gene.ref,Gene,Description)

paths.summary[paths.summary$num_nodes < 20 |paths.summary$num_edges <20,]
paths.summary[paths.summary$eigen >10,]
paths.summary[paths.summary$num_nodes > 1000 |paths.summary$num_edges >4000,]

all.path.names     <- pathway.titles
all.hsap.essential <- data_frame()

saveRDS(paths.summary,"human_data/paths_summary.RDS")
saveRDS(all.path.names,"human_data/all_path_names.RDS")
saveRDS(graphs.hsap,"human_data/graphs_hsap.RDS")
saveRDS(hsap.gene.ref,"human_data/hsap_gene_ref.RDS")
