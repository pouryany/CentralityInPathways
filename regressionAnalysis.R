library(dplyr)
library(broom)
library(ggplot2)
library(xtable)

gene.essential <- readRDS("gene_essentials.rds")


## The regression analysis

a <- gene.essential %>% group_by(katz.ssc.quant, Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer") %>%
    mutate(., Centrality = "Katz-Source-Sink",quant = katz.ssc.quant)

b <- gene.essential %>% group_by(deg.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer") %>%
    mutate(Centrality = "Degree",quant = deg.quant )

c <- gene.essential %>% group_by(beet.source.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer") %>%
    mutate(Centrality = "Betweenness" ,quant = beet.source.quant)

d <- gene.essential %>% group_by(pgr.source.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer") %>%
    mutate(Centrality = "PageRank",quant = pgr.source.quant)

e <- gene.essential %>% group_by(katz.source.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "Katz",quant = katz.source.quant)

f <- gene.essential %>% group_by(pgr.ssc.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "PageRank-Source-Sink",quant = pgr.ssc.quant )

g <- gene.essential %>% group_by(pgr.und.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "PageRank-Undirected",quant = pgr.und.quant )

h <- gene.essential %>% group_by(lap.source.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "Lap-SSC",quant = lap.source.quant )

i <- gene.essential %>% group_by(beet.und.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "Beetween",quant = beet.und.quant )

total <- rbind(a,b,d,e,f,g,h,i)



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

