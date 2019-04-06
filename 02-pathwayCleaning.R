library(dplyr)
library(ggplot2)
all.homo.essential <- readRDS("pathwayCentralities.rds")

gene.essential <-  all.homo.essential
gene.essential %>% distinct(.,pathway.name)
gene.essential <-  gene.essential %>%  replace_na(list(Description = "Normal"))

# Assigning quantile scores to centrality models
gene.essential %<>% filter(.,total.node < 1000, total.node >20, total.edge > 20,
                           total.edge < 4000) %>%
    mutate(., katz.ssc.quant =    round((katz.ssc.rank/total.node)*100),
           katz.und.quant =    round((katz.und.rank/total.node)*100),
           katz.sink.quant =   round((katz.sink.rank/total.node)*100),
           katz.source.quant = round((katz.source.rank/total.node)*100),
           deg.quant =         round((degree.rank/total.node)*100),
           beet.source.quant = round((beet.source.rank/total.node)*100),
           beet.und.quant =    round((beet.und.rank/total.node)*100),
           beet.sink.quant =   round((beet.sink.rank/total.node)*100),
           beet.ssc.quant =    round((beet.ssc.rank/total.node)*100),
           pgr.source.quant =  round((pgr.source.rank/total.node)*100),
           pgr.sink.quant =    round((pgr.sink.rank/total.node)*100),
           pgr.ssc.quant =     round((pgr.ssc.rank/total.node)*100),
           pgr.und.quant =     round((pgr.und.rank/total.node)*100),
           lap.ssc.quant =     round((lap.ssc.rank/total.node)*100),
           lap.source.quant =  round((lap.source.rank/total.node)*100),
           lap.sink.quant =    round((lap.sink.rank/total.node)*100),
           lap.und.quant =     round((lap.und.rank/total.node)*100))


### A total of 219 pathways passed the criteria
gene.essential %>% distinct(.,pathway.name)



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



gene.essential <- gene.essential[ !(gene.essential$pathway.name %in% unlist(no.cancers))  ,]


saveRDS(gene.essential, file = "gene_essentials.rds")

