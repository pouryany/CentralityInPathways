library(graphite)
library(stringr)
source('https://bioconductor.org/biocLite.R')
biocLite('rentrez')
library(rentrez)
library(Rgraphviz)
library("org.Hs.eg.db")





if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("org.Hs.egUNIGENE.db", version = "3.8")



select(org.Hs.eg.db, temp, "SYMBOL", "UNIPROT")

?select

rm(list = ls())

pathlist <- pathways("hsapiens", "reactome")



pathGraphs <- sapply(pathlist, pathwayGraph)


unigene <- toTable(org.Hs.egUNIGENE)
# extract 100 random unigene entries
id.type  <- "unigene"
res <- getSYMBOL(temp,data = "org.Hs.egUNIGENE.db") 


temp <- nodes(pathGraphs[[8]])
temp <- str_sub(temp, start = 9)
temp <- getSYMBOL(temp, data='org.Hs.eg.db')
nodes(X) <- unname(temp)

getSYMBOL("16197", data='org.Hs.egUNIGENE.db')


graphSizes <- sapply(pathGraphs, function(X)({
    length(nodes(X))
}))


selvec <- graphSizes > 20 & graphSizes < 1000

pathGraphs <- pathGraphs[selvec]

length(pathGraphs)

testPathway <- readRDS("testpathway.RDS")


plot(testPathway)
graphs.homo    <-  sapply(pathGraphs, function(X) ({
    temp <- nodes(X)
    temp <- str_sub(temp, start = 10)
    temp <- getSYMBOL(temp, data='org.Hs.eg')
    if(sum(is.na(temp)) > 0){
        no.names <- paste0("NoName",1:sum(is.na(temp)))
        temp[is.na(temp)] <- no.names
    }
    
    nodes(X) <- unname(temp)
    return(X)
}))

i <-0
for (X in pathGraphs){
    temp <- nodes(X)
    temp <- str_sub(temp, start = 10)
    temp <- getSYMBOL(temp, data='org.Hs.eg')
    nodes(X) <- unname(temp)
    print(i)
    i <- i+1
}


?getSYMBOL()

 nodes(pathGraphs[[3]])
 
 
 
 
 
 
 
 
 ### Graphite sucks, lets try CePA
 ## Also check this article 
 ## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4009765/
 
 BiocManager::install("CePa", version = "3.8")
 rm(list = ls())
 library(CePa)
 
 data("PID.db")
 
reactome <- PID.db$Reactome
biocarta <- PID.db$BioCarta
biocarta$pathList
class(biocarta)
get.cepa(biocarta)

plot(PID.db$NCI)
