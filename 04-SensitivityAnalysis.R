# We are investigating how much error in the edge annotations affects the 
# rankings from a centrality model. So let's investigate. 

# Let's take one pathway. 

library(GGally)
library(network)
library(sna)
library(ggplot2)
library(RBGL)
library(KEGGgraph)
library(org.Hs.eg.db)
library(annotate)
library(Hmisc)
library(CADIA)

pathways.collection.names
## Getting a sample pathway (ErbB) from CADIA prefiltered dataset.


## ErbB is "04012.xml  Ras is "04014.xml
apop.Graph   <- pathways.collection[["04012.xml"]]
node.names   <- nodes(apop.Graph)
node.ENTREZ  <- KEGGgraph::translateKEGGID2GeneID(node.names)
node.genes   <- getSYMBOL(node.ENTREZ,data='org.Hs.eg')
apop.mat     <- as(apop.Graph,"matrix")



## Now we randomly add an non existing edge and remove it to see how much 
## The centrality chagnes. Let's see how do we do it with the matrix


### Interesting observation: As I increase alpha something happens:
### therese is a cubic component, then it becomes linear and then it becomes
### constant.  The higher alpha, these events happen faster. 
### https://www.cs.cmu.edu/~avrim/598/chap4only.pdf This link might help.



# Enter Source-Sink Centrality
#apop.mat <- matrix(0,nrow(apop.mat),nrow(apop.mat))

ssc.mat  <- CADIA::newpath.centrality(apop.mat,alpha = 0.2, beta = 1)
ssc.vec  <- rowSums(ssc.mat)
ssc.rank <- rank(ssc.vec, ties.method = "min")



#
pert.range  <- 300
pert.rpt    <- 50

dif.array <- matrix(NA,pert.rpt,pert.range)
zero.inds <- which(apop.mat ==0)
for (j in 1:pert.range){
    for (i in 1:pert.rpt){
        set.seed((10000*j)+i)
        temp.inds <- sample(zero.inds,j)
        apop.mat[temp.inds] <- 1
        cent.rand.mat  <- CADIA::newpath.centrality(apop.mat,alpha = 0.2,
                                                    beta = 1)
        cent.rand.vec  <- rowSums(cent.rand.mat)
        rand.rank      <- rank(cent.rand.vec, ties.method = "min")
        dif.array[i,j]   <- sqrt(sum((ssc.rank - rand.rank) **2)) 
        
        apop.mat[temp.inds] <- 0
    } 
}



### This observation is fantastic!!! Difference is cubic with number of edges
plot(1:pert.range,(apply(dif.array, 2, mean)))


dif.array[,70]





katz.mat  <- CADIA::newpath.centrality(apop.mat,alpha = 0.1, beta = 0)
katz.vec  <- rowSums(katz.mat)
katz.rank <- rank(katz.vec, ties.method = "min")


dif.array.katz <- matrix(NA,100,200)
zero.inds <- which(apop.mat ==0)
for (j in 1:200){
    for (i in 1:100){
        set.seed((10000*j)+i)
        temp.inds <- sample(zero.inds,j)
        apop.mat[temp.inds] <- 1
        cent.rand.mat  <- CADIA::newpath.centrality(apop.mat,alpha = 0.1,
                                                    beta = 0)
        cent.rand.vec  <- rowSums(cent.rand.mat)
        rand.rank      <- rank(cent.rand.vec, ties.method = "min")
        dif.array.katz[i,j]   <- sqrt(sum((katz.rank - rand.rank) **2)) 
        
        apop.mat[temp.inds] <- 0
    } 
}

plot(1:200,(apply(dif.array.katz, 2, mean))**2.9)



### On with PageRank Source-Sink 



#PageRank Source-Sink Centrality
pgr.ssc <- function(someMat){
    igraph.obj   <- igraph::graph_from_adjacency_matrix(someMat,mode = "directed")
    igraph.obj2  <- igraph::graph_from_adjacency_matrix(t(someMat),mode = "directed")
    
    ### Source Directed pageRank
    pgr.source.vec   <- igraph::page.rank(igraph.obj,damping = 0.9)
    pgr.source.vec   <- pgr.source.vec$vector
    
    ###Source-Sink pageRank
    pgr.sink.vec     <- igraph::page.rank(igraph.obj2,damping = 0.9)
    pgr.sink.vec     <- pgr.sink.vec$vector
    pgr.ssc.vec      <- pgr.source.vec + pgr.sink.vec
    
    pgr.ssc.rank <- rank(pgr.ssc.vec, ties.method = "min")
    return(pgr.ssc.rank)
}


apop.pgr.rank <- pgr.ssc(apop.mat)

dif.array.pgssc <- matrix(NA,500,200)
zero.inds <- which(apop.mat ==0)
for (j in 1:200){
    for (i in 1:500){
        set.seed((10000*j)+i)
        temp.inds <- sample(zero.inds,j)
        apop.mat[temp.inds] <- 1
        cent.rand.vec  <- pgr.ssc(apop.mat)
        rand.rank      <- rank(cent.rand.vec, ties.method = "min")
        dif.array.pgssc[i,j]   <- sqrt(sum((apop.pgr.rank - rand.rank) **2)) 
        
        apop.mat[temp.inds] <- 0
    } 
}

plot(1:200,(apply(dif.array.pgssc, 2, mean))**5)




