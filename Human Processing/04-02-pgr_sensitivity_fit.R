rm(list = ls())
# We first calculate the PageRank centrality variations for a grid of alpha
# We then do the regression analysis.


library(dplyr)
library(ggplot2)
library(broom)
library(xtable)
library(reshape2)

all.hsap.essential <- readRDS("human_data/pgr_sensitivity_analysis.RDS")

gene.essential <-  all.hsap.essential
gene.essential %>% distinct(.,pathway.name)
gene.essential <-  gene.essential %>%  tidyr::replace_na(list(Description = "Normal"))

# Assigning quantile scores to centrality models
gene.essential %<>% filter(.,total.node < 1000, total.node >20, total.edge > 20,
                           total.edge < 4000) %>%
                   mutate(.,
                           pgr.source.quant =  round((pgr.source.rank/total.node)*100),
                           pgr.sink.quant =    round((pgr.sink.rank/total.node)*100),
                           pgr.ssc.quant =     round((pgr.ssc.rank/total.node)*100),
                           pgr.und.quant =     round((pgr.und.rank/total.node)*100))


### A total of 157 pathways passed the criteria
gene.essential %>% distinct(.,pathway.name)



gene.essential$Description <- factor(gene.essential$Description)
gene.essential$Description <- relevel(gene.essential$Description,ref = "Normal")



feature.list <- grep("quant",names(gene.essential))
new.vars     <- c("pathway.name","node.genes","Description","pgr_damp",
                  names(gene.essential)[feature.list])

cents.table  <- gene.essential[,new.vars]
cents.table  <- reshape2::melt(cents.table, id.vars = c("pathway.name",
                                                        "node.genes","Description","pgr_damp"))
names(cents.table)[5] <- "Centrality"


#cents.table <- cents.table[cents.table$pgr_damp == 0.86,]
tots <- cents.table %>% dplyr::group_by(.,pgr_damp,Centrality,value,Description)  %>%
    dplyr::summarise(n = dplyr::n()) %>% dplyr::mutate(freq = n/ sum(n)) %>%
    dplyr::filter(., Description== "Cancer") 


overall.cors <- tots %>% group_by(.,pgr_damp,Centrality) %>%
    do(fit = lm(freq ~ as.numeric(value), data=.)) %>% 
    glance(fit) 






total <- overall.cors

library(RColorBrewer)
colores <- brewer.pal(4,"Set2")  
ggplot(total, aes(y = adj.r.squared, x= pgr_damp)) + 
    geom_line(aes(color = Centrality), size = 2) +
    scale_color_manual( labels = c("Source","Sink","SSC","Undirected"),
                        values = c("pgr.source.quant" = colores[1],
                                   "pgr.sink.quant" = colores[2],
                                   "pgr.ssc.quant" = colores[3],
                                   "pgr.und.quant" = colores[4])) +
    labs(x = "PageRank dampening factor (alpha)", 
         y = "Adjusted R-Squared of the linear model")+
    theme_bw()+
    theme(strip.text = element_text(face="bold", size=20),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.title=element_text(face = "bold", size = 20),
          axis.text.y=element_text(size = 12),
          axis.text.x=element_text(size = 12),
          axis.ticks.y=element_blank())






overall.cors <- tots %>% group_by(Centrality) %>%
    do(fit = lm(freq ~ as.numeric(value), data=.)) %>% 
    glance(fit) %>% filter(., grepl("cls",Centrality))



overall.cors     <- overall.cors[,c(1,2,3,6)]
adjusted.pval    <- p.adjust(unlist(overall.cors[,4]),method = "fdr")
adjusted.pval    <- formatC(adjusted.pval, format = "e", digits = 2)
overall.cors[,4] <- formatC(unlist(pull(overall.cors[,4])), format = "e",
                            digits = 2)


print(xtable(overall.cors), include.rownames = F)

text.vals <- paste("Adj r-squared:",formatC(pull(overall.cors[,2]),digits = 2),
                   "\n p-value:", pull(overall.cors[,4]), "\n Adj p-value:",
                   adjusted.pval)

overall.cors$label <- text.vals


total <- tots %>% filter(.,grepl("cls",Centrality))

ggplot(total, aes(y = 100*freq, x= value)) + geom_point()+
    geom_smooth(method= "lm") + #geom_smooth(method= "loess", color="green" , fill = "red") +
    facet_wrap(~Centrality ,ncol = 2) +theme_bw()+
    labs(x = "Normalized Centrality Score (Eq. 20)", y = "% of Genes that are Cancer-related (Eq. 21)")+
    theme(strip.text = element_text(face="bold", size=20),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 30),
          legend.text = element_text(size = 9),
          legend.title=element_text(face = "bold", size = 9),
          axis.text.y=element_text(size = 12),
          axis.text.x=element_text(size = 12),
          axis.ticks.y=element_blank()) +
    geom_text(data=overall.cors, x = 50, y =70,
              aes(x,y,label=label),size = 6, inherit.aes=FALSE)
ggsave("images/Human_Regression_Closeness.pdf",
       width = 10, height = 10, units = c("in"))


