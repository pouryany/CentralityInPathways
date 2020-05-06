rm(list = ls())
library(dplyr)
library(broom)
library(ggplot2)
library(xtable)
library(reshape2)


gene.essential <- readRDS("human_data/gene_essentials.rds")
gene.essential$Description <- factor(gene.essential$Description)
gene.essential$Description <- relevel(gene.essential$Description,ref = "Normal")

path.items <- unique(gene.essential$pathway.name)
saveRDS(path.items,"human_data/pathItems.RDS")


feature.list <- grep("quant",names(gene.essential))
new.vars     <- c("pathway.name","node.genes","Description",
                   names(gene.essential)[feature.list])
new.vars     <- new.vars[-c(5)]
cents.table  <- gene.essential[,new.vars]
cents.table  <- reshape2::melt(cents.table, id.vars = c("pathway.name",
                                                  "node.genes","Description"))
names(cents.table)[4] <- "Centrality"


tots <- cents.table %>%  dplyr::group_by(.,Centrality,value,Description)  %>%
                         dplyr::summarise(n = dplyr::n()) %>% dplyr::mutate(freq = n/ sum(n)) %>%
                         dplyr::filter(., Description== "Cancer") 


overall.cors <- tots %>% group_by(Centrality) %>%
                do(fit = lm(freq ~ as.numeric(value), data=.)) %>% 
                glance(fit) %>% filter(., !grepl("cls",Centrality))



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


total <- tots %>% filter(.,!grepl("cls",Centrality))

ggplot(total, aes(y = 100*freq, x= value)) + geom_point()+
  geom_smooth(method= "lm") + #geom_smooth(method= "loess", color="green" , fill = "red") +
  facet_wrap(~Centrality ,ncol = 3) +theme_bw()+
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
ggsave("images/Human_Regression.pdf",
       width = 10, height = 18, units = c("in"))





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


total  <- tots %>% filter(.,grepl("cls",Centrality))
x.name <- levels(total$Centrality)

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


if(FALSE){
  ## The regression analysis
  
  a1 <- gene.essential %>% group_by(katz.ssc.quant, Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer") %>%
    mutate(., Centrality = "Katz-SSC",quant = katz.ssc.quant)
  
  a2 <- gene.essential %>% group_by(katz.source.quant, Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer") %>%
    mutate(., Centrality = "Katz-Source",quant = katz.source.quant)
  
  a3 <- gene.essential %>% group_by(katz.sink.quant, Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer") %>%
    mutate(., Centrality = "Katz-Sink",quant = katz.sink.quant)
  
  b <- gene.essential %>% group_by(deg.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer") %>%
    mutate(Centrality = "Degree",quant = deg.quant )
  
  d1 <- gene.essential %>% group_by(pgr.source.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer") %>%
    mutate(Centrality = "PageRank-Source",quant = pgr.source.quant)
  
  d2 <- gene.essential %>% group_by(pgr.sink.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer") %>%
    mutate(Centrality = "PageRank-Sink",quant = pgr.sink.quant)
  
  d3 <- gene.essential %>% group_by(pgr.ssc.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer") %>%
    mutate(Centrality = "PageRank-SSC",quant = pgr.ssc.quant)
  
  
  d4 <- gene.essential %>% group_by(pgr.und.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "PageRank-Und",quant = pgr.und.quant )
  
  h1 <- gene.essential %>% group_by(lap.source.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "Lap-Source",quant = lap.source.quant )
  
  h2 <- gene.essential %>% group_by(lap.sink.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "Lap-Sink",quant = lap.sink.quant )
  
  h3 <- gene.essential %>% group_by(lap.ssc.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "Lap-SSC",quant = lap.ssc.quant )
  
  
  # Closeness centrality 
  g1 <- gene.essential %>% group_by(cls.source.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "Closeness-SSC",quant = cls.source.quant )
  
  g2 <- gene.essential %>% group_by(cls.sink.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "Closeness-SSC",quant = cls.sink.quant )
  
  g3 <- gene.essential %>% group_by(cls.und.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "Closeness-SSC",quant = cls.und.quant )
  
  g4 <- gene.essential %>% group_by(cls.ssc.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "Closeness-SSC",quant = cls.ssc.quant )
  
  
  total <- rbind(a1,a2,a3,b,d1,d2,d3,d4,h1,h2,h3)
  
  tot2  <- total[total$Centrality %in% c("PageRank-Und","PageRank-SSC"),]
  #rm(a1,a2,a3,b,c1,c2,c3,c4,d1,d2,d3,d4,h1,h2,h3,h4)
  
  
  # Reporting
  
  overall.cors <- total %>% group_by(Centrality) %>%
    do(fit = lm(freq ~ as.numeric(quant), data=.)) %>%  glance(fit)
  
  
  
  overall.cors  <- overall.cors[,c(1,2,3,6)]
  adjusted.pval <- p.adjust(unlist(overall.cors[,4]),method = "fdr")
  adjusted.pval <- formatC(adjusted.pval, format = "e", digits = 2)
  overall.cors[,4] <- formatC(unlist(pull(overall.cors[,4])), format = "e", digits = 2)
  print(xtable(overall.cors), include.rownames = F)
  
  text.vals <- paste("Adj r-squared:",formatC(pull(overall.cors[,2]),digits = 2),
                     "\n p-value:", pull(overall.cors[,4]), "\n Adj p-value:",
                     adjusted.pval)
  
  overall.cors$label <- text.vals
  
  ggplot(total, aes(y = 100*freq, x= quant)) + geom_point()+
    geom_smooth(method= "lm") + #geom_smooth(method= "loess", color="green" , fill = "red") +
    facet_wrap(~Centrality ,ncol = 3) +theme_bw()+
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
  
  ggsave("images/Human_Regression.pdf",
         width = 10, height = 18, units = c("in"))
  
  
  # Calculating precentage
  
  overall.regs <- total %>% group_by(Centrality) %>%
    do(fit = lm((100*freq) ~ as.numeric(quant), data=.)) %>%  tidy(fit)
  
  
  overall.regs$term[overall.regs$term == "as.numeric(quant)"] <- "Coefficient"
  
  overall.regs[,c(3,4,5,6)] <- sapply(overall.regs[,c(3,4,5,6)], FUN = formatC,
                                      format = "e", digits = 2
  )
  
  
  print(xtable(overall.regs), include.rownames = F, digits = 2,
        display=c("s","s", "s", "s","g"),
        math.style.exponents = T)
  
  
  tot2  <- total[total$Centrality %in% c("PageRank-Und","PageRank-SSC"),]
  
  m.interaction <- lm((100*freq) ~ as.numeric(quant)/Centrality, data=tot2)
  summary(m.interaction)
  ### Functionalize the above code at some point. Use the chunks below
  
  
  
  
  
  
  
  
  
  # Closeness centrality 
  g1 <- gene.essential %>% group_by(cls.source.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "Closeness-Source",quant = cls.source.quant )
  
  g2 <- gene.essential %>% group_by(cls.sink.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "Closeness-Sink",quant = cls.sink.quant )
  
  g3 <- gene.essential %>% group_by(cls.und.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "Closeness-Und",quant = cls.und.quant )
  
  g4 <- gene.essential %>% group_by(cls.ssc.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Cancer")  %>%
    mutate(Centrality = "Closeness-SSC",quant = cls.ssc.quant )
  
  total <- rbind(g1,g2,g3,g4)
  
  #rm(a1,a2,a3,b,c1,c2,c3,c4,d1,d2,d3,d4,h1,h2,h3,h4)
  
  
  # Reporting
  
  overall.cors <- total %>% group_by(Centrality) %>%
    do(fit = lm(freq ~ as.numeric(quant), data=.)) %>%  glance(fit)
  
  
  
  overall.cors  <- overall.cors[,c(1,2,3,6)]
  adjusted.pval <- p.adjust(unlist(overall.cors[,4]),method = "fdr")
  adjusted.pval <- formatC(adjusted.pval, format = "e", digits = 2)
  overall.cors[,4] <- formatC(unlist(pull(overall.cors[,4])), format = "e", digits = 2)
  print(xtable(overall.cors), include.rownames = F)
  
  text.vals <- paste("Adj r-squared:",formatC(pull(overall.cors[,2]),digits = 2),
                     "\n p-value:", pull(overall.cors[,4]), "\n Adj p-value:",
                     adjusted.pval)
  
  overall.cors$label <- text.vals
  
  ggplot(total, aes(y = 100*freq, x= quant)) + geom_point()+
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
  
  ggsave("images/Human_Regression_closeness.pdf",
         width = 10, height = 12, units = c("in"))
  
  
  # Calculating precentage
  
  overall.regs <- total %>% group_by(Centrality) %>%
    do(fit = lm((100*freq) ~ as.numeric(quant), data=.)) %>%  tidy(fit)
  
  
  overall.regs$term[overall.regs$term == "as.numeric(quant)"] <- "Coefficient"
  
  overall.regs[,c(3,4,5,6)] <- sapply(overall.regs[,c(3,4,5,6)], FUN = formatC,
                                      format = "e", digits = 2
  )
  
  
  print(xtable(overall.regs), include.rownames = F, digits = 2,
        display=c("s","s", "s", "s","g"),
        math.style.exponents = T)
  
  
  tot2  <- total[total$Centrality %in% c("PageRank-Und","PageRank-SSC"),]
  
  m.interaction <- lm((100*freq) ~ as.numeric(quant)/Centrality, data=tot2)
  summary(m.interaction)
  
  # 
  # var  <- "katz.ssc.quant"
  # cent <- "Katz_Source_Sink"
  # Desc <- "Description"
  # 
  # quant.sorter <- function(data,var, Desc,cent){
  # 
  #     varval1 <- lazyeval::interp(~quant, quant = as.name(var))
  #     varval2 <- lazyeval::interp(~Centrality, Centrality = cent)
  #     sorted  <- gene.essential %>% group_by_(var, as.name(Desc))  %>%
  #         summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
  #         filter_(., ~Description == "Cancer") %>%
  #         mutate_(.dots = setNames(list(varval1,varval2),
  #                                  c("quant", "Centrality")))
  # 
  #     return(sorted)
  # }
  # 
  # 
  # model.vars <-  c("katz.ssc.quant", "katz.source.quant","katz.sink.quant",
  #                  "deg.quant","beet.source.quant","beet.sink.quant",
  #                  "beet.ssc.quant","beet.und.quant","pgr.source.quant",
  #                  "pgr.sink.quant","pgr.ssc.quant","pgr.und.quant",
  #                  "lap.source.quant","lap.sink.quant","lap.ssc.quant")
  # 
  # model.names <-  c("katz_ssc", "katz_source","katz_sink",
  #                   "deg","beet_source","beet_sink",
  #                   "beet_ssc","beet_und","pgr_source",
  #                   "pgr_sink","pgr_ssc","pgr_und",
  #                   "lap_source","lap_sink","lap_ssc")
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  
}

