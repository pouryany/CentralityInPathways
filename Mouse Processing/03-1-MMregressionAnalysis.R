rm(list = ls())
library(dplyr)
library(broom)
library(ggplot2)
library(xtable)
library(reshape2)

gene.essential <- readRDS("mouse_data/MMgene_essentials.rds")

gene.essential$Description <- factor(gene.essential$Description)
gene.essential$Description <- relevel(gene.essential$Description,ref = "Normal")


feature.list <- grep("quant",names(gene.essential))
new.vars     <- c("pathway.name","node.genes","Description",
                  names(gene.essential)[feature.list])
new.vars     <- new.vars[-c(5)]
cents.table  <- gene.essential[,new.vars]
cents.table  <- reshape2::melt(cents.table, id.vars = c("pathway.name",
                                                        "node.genes","Description"))
names(cents.table)[4] <- "Centrality"


tots <- cents.table %>%  group_by(.,Centrality,value,Description)  %>%
  summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
  filter(., Description== "Lethal") 


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

overall.cors$labels <- text.vals


total <- tots %>% filter(.,!grepl("cls",Centrality))

ggplot(total, aes(y = 100*freq, x= value)) + geom_point()+
  geom_smooth(method= "lm") + #geom_smooth(method= "loess", color="green" , fill = "red") +
  facet_wrap(~Centrality ,ncol = 3) +theme_bw()+
  labs(x = "Normalized Centrality Score (Eq. 20)", y = "% of Genes that are Lethal (Eq. 21)")+
  theme(strip.text = element_text(face="bold", size=20),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 30),
        legend.text = element_text(size = 9),
        legend.title=element_text(face = "bold", size = 9),
        axis.text.y=element_text(size = 12),
        axis.text.x=element_text(size = 12),
        axis.ticks.y=element_blank()) +
  geom_text(data=overall.cors, x = 50, y =70,
            aes(x,y,label=labels),size = 6, inherit.aes=FALSE)
ggsave("images/Mouse_Regression.pdf",
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


total <- tots %>% filter(.,grepl("cls",Centrality))

ggplot(total, aes(y = 100*freq, x= value)) + geom_point()+
  geom_smooth(method= "lm") + #geom_smooth(method= "loess", color="green" , fill = "red") +
  facet_wrap(~Centrality ,ncol = 2) +theme_bw()+
  labs(x = "Normalized Centrality Score (Eq. 20)", y = "% of Genes that are Lethal (Eq. 21)")+
  theme(strip.text = element_text(face="bold", size=20),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 30),
        legend.text = element_text(size = 9),
        legend.title=element_text(face = "bold", size = 9),
        axis.text.y=element_text(size = 12),
        axis.text.x=element_text(size = 12),
        axis.ticks.y=element_blank()) +
  geom_text(data=overall.cors, x = 50, y =50,
            aes(x,y,label=label),size = 6, inherit.aes=FALSE)
ggsave("images/Mouse_Regression_Closeness.pdf",
       width = 10, height = 10, units = c("in"))











## The regression analysis
if(FALSE){
a1 <- gene.essential %>% group_by(katz.ssc.quant, Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Lethal") %>%
    mutate(., Centrality = "Katz-Source-Sink",quant = katz.ssc.quant)

a2 <- gene.essential %>% group_by(katz.source.quant, Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Lethal") %>%
    mutate(., Centrality = "Katz-Source",quant = katz.source.quant)

a3 <- gene.essential %>% group_by(katz.sink.quant, Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Lethal") %>%
    mutate(., Centrality = "Katz-Sink",quant = katz.sink.quant)

b <- gene.essential %>% group_by(deg.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Lethal") %>%
    mutate(Centrality = "Degree",quant = deg.quant )



d1 <- gene.essential %>% group_by(pgr.source.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Lethal") %>%
    mutate(Centrality = "PageRank-Source",quant = pgr.source.quant)

d2 <- gene.essential %>% group_by(pgr.sink.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Lethal") %>%
    mutate(Centrality = "PageRank-Sink",quant = pgr.sink.quant)

d3 <- gene.essential %>% group_by(pgr.ssc.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Lethal") %>%
    mutate(Centrality = "PageRank-SSC",quant = pgr.ssc.quant)


d4 <- gene.essential %>% group_by(pgr.und.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Lethal")  %>%
    mutate(Centrality = "PageRank-Und",quant = pgr.und.quant )

h1 <- gene.essential %>% group_by(lap.source.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Lethal")  %>%
    mutate(Centrality = "Lap-Source",quant = lap.source.quant )

h2 <- gene.essential %>% group_by(lap.sink.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Lethal")  %>%
    mutate(Centrality = "Lap-Sink",quant = lap.sink.quant )

h3 <- gene.essential %>% group_by(lap.ssc.quant,Description)  %>%
    summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Lethal")  %>%
    mutate(Centrality = "Lap-SSC",quant = lap.ssc.quant )



total <- rbind(a1,a2,a3,b,d1,d2,d3,d4,h1,h2,h3)

#rm(a1,a2,a3,b,c1,c2,c3,c4,d1,d2,d3,d4,h1,h2,h3,h4)


# Reporting


overall.cors <- total %>% group_by(Centrality) %>%
    do(fit = lm(freq ~ as.numeric(quant), data=.)) %>%  glance(fit)



overall.cors  <- overall.cors[,c(1,2,3,6)]
adjusted.pval <- p.adjust(unlist(overall.cors[,4]),method = "fdr")
adjusted.pval <- formatC(adjusted.pval, format = "e", digits = 2)
overall.cors[,4] <- formatC(unlist(pull(overall.cors[,4])), format = "e", digits = 2)
print(xtable(overall.cors), include.rownames = F)

text.vals <- paste("Adjusted r-squared:",formatC(pull(overall.cors[,2]),digits = 2),
                   "\n p-value:", pull(overall.cors[,4]), "\n Adj p-value:",
                   adjusted.pval)

overall.cors$label <- text.vals


ggplot(total, aes(y = 100 * freq, x= quant)) + geom_point()+
    geom_smooth(method= "lm") + #geom_smooth(method= "loess", color="green" , fill = "red") +
    facet_wrap(~Centrality ,ncol = 3) +theme_bw()+
    labs(x = "Normalized Centrality Score (Eq. 20)", y = "% of Genes that are Lethal (Eq. 21)")+
    theme(strip.text = element_text(face="bold", size=20),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 30),
          legend.text = element_text(size = 9),
          legend.title=element_text(face = "bold", size = 9),
          axis.text.y=element_text(size = 12),
          axis.text.x=element_text(size = 12),
          axis.ticks.y=element_blank()) +
    geom_text(data=overall.cors, x = 50, y =30,
              aes(x,y,label=label),size = 6, inherit.aes=FALSE)

ggsave("images/Mouse_Regression.pdf",
       width = 10, height = 18, units = c("in"))


#reporting percentage
overall.regs <- total %>% group_by(Centrality) %>%
    do(fit = lm((100*freq) ~ as.numeric(quant), data=.)) %>%  tidy(fit)


overall.regs$term[overall.regs$term == "as.numeric(quant)"] <- "Coefficient"

overall.regs[,c(3,4,5,6)] <- sapply(overall.regs[,c(3,4,5,6)], FUN = formatC,
                                    format = "e", digits = 2
                                    )


print(xtable(overall.regs), include.rownames = F, digits = 2,
      display=c("s","s", "s", "s","g"),
      math.style.exponents = T)


# Closeness centrality 
g1 <- gene.essential %>% group_by(cls.source.quant,Description)  %>%
  summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
  filter(., Description== "Lethal")  %>%
  mutate(Centrality = "Closeness-Source",quant = cls.source.quant )

g2 <- gene.essential %>% group_by(cls.sink.quant,Description)  %>%
  summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
  filter(., Description== "Lethal")  %>%
  mutate(Centrality = "Closeness-Sink",quant = cls.sink.quant )

g3 <- gene.essential %>% group_by(cls.und.quant,Description)  %>%
  summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
  filter(., Description== "Lethal")  %>%
  mutate(Centrality = "Closeness-Und",quant = cls.und.quant )

g4 <- gene.essential %>% group_by(cls.ssc.quant,Description)  %>%
  summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
  filter(., Description== "Lethal")  %>%
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

ggsave("images/Mouse_Regression_closeness.pdf",
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

}

