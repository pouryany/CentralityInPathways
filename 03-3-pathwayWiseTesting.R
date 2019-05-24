###
###  Pathway wise testing
###  Gene essential is cleaned and without cancer-less pathways
rm(list = ls())
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(xtable)

gene.essential <- readRDS("gene_essentials.rds")
# Pathway-wise checking if the rank of the cancer nodes are different based on centrality
### Getting rid of pathways with no cancer-genes


p.vals.process <- gene.essential
distinct(p.vals.process,pathway.name)
# p.vals.process <- p.vals.process %>% mutate(.,pgr.log.vec2 = log(pgr.ssc.vec),
#                                             cent.log.vec = log(katz.ssc.vec),
#                                             pgr.log.vec = log(pgr.source.vec),
#                                             pgr.log.vec3= log(pgr.und.vec))


# You can change the t.test to wilcox.test in the below formula to compute
# Nonparametric statistics. In that case, log values would not differ.
# For publication remember to add the other possible alternative hypothesis.
p.vals.process <- gene.essential %>%
    gather(., key = "Centrality", value = "cent_value",
           pgr.ssc.vec,pgr.source.vec,pgr.sink.vec,pgr.und.vec,
           degree.norm,
           katz.ssc.vec,katz.source.norm,katz.sink.norm,
          lap.ssc.vec,lap.sink.vec,lap.source.vec)%>%
    group_by(pathway.name,Centrality) %>%
    do(pval =tidy(wilcox.test( cent_value~ Description,
                          alternative = "greater",paired = F,exact=FALSE, data = .))) %>%
    unnest()  %>%
    group_by(Centrality) %>%
    mutate(., fdr = p.adjust(p.value)) %>%
    filter(., fdr < 0.05)


# This plot is intended to investigate false positives. Include later
# You have to take out the fdr filter from the above if you want it work

ggplot(p.vals.process, aes(x= p.value)) +
    geom_density(adjust = 1/4) + #geom_smooth(method= "loess", color="green" , fill = "red") +
    facet_wrap(~Centrality ,ncol = 4) +theme_bw()+
    theme(strip.text = element_text(face="bold", size=20),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 30),
          legend.text = element_text(size = 9),
          legend.title=element_text(face = "bold", size = 9),
          axis.text.y=element_text(size = 12),
          axis.text.x=element_text(size = 12),
          axis.ticks.y=element_blank())


# Creating confusion matrix
mat1 <- table(p.vals.process$Centrality,p.vals.process$pathway.name)
mat2 <- table(p.vals.process$pathway.name,p.vals.process$Centrality)


confusion <- as.matrix(mat1) %*% as.matrix(mat2)

colnames(confusion) <- c("Degree", "Katz-Sink", "Katz-Source","Katz-SSC",
                        "Lap-Sink","Lap-Source","Lap-SSC",
                        "Pgr-Sink","Pgr-Source", "Pgr-SSC", "Pgr-Und")
rownames(confusion) <- colnames(confusion)

# The following line only for nonparametrix you have to remove the log values
# The Wilcox test will not make any difference with the log values
#confusion2 <- confusion[c("SS-Katz","Degree", "Katz", "PageRank",
#                          "SS-PageRank", "Und. PageRank"),
#                        c("SS-Katz","Degree", "Katz", "PageRank",
#                          "SS-PageRank", "Und. PageRank")]

#confusion2[lower.tri(confusion2,diag = F)] <- ""
confusion[lower.tri(confusion,diag = F)] <- ""
xtable(confusion)

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
gene.essential %>% distinct(., node.genes)
gene.essential %>% filter(.,Description == "Cancer") %>% distinct(., node.genes)

