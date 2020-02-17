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


library(grid)
# Create a text
grobA <- grobTree(textGrob("A", x=0.05,  y=0.95, hjust=0,
                          gp=gpar( fontsize=30, fontface="bold")))

grobB <- grobTree(textGrob("B", x=0.05,  y=0.95, hjust=0,
                           gp=gpar( fontsize=30, fontface="bold")))
# Plot


library(RColorBrewer)
colores <- brewer.pal(4,"Set2")  
p1 <- ggplot(total, aes(y = adj.r.squared, x= pgr_damp)) + 
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
              axis.text.y=element_text(size = 15),
              axis.text.x=element_text(size = 15),
              axis.ticks.y=element_blank())+ 
    annotation_custom(grobA)

    






max.r <- which.max(overall.cors$adj.r.squared)
max.r <- overall.cors[max.r,]$pgr_damp











overall.cors <- tots %>% group_by(.,pgr_damp,Centrality) %>%
    do(fit = lm(freq ~ as.numeric(value), data=.)) %>% 
    glance(fit) %>% filter(.,pgr_damp == max.r)


psych::r.test(100,0.808,0.692, twotailed = F)

overall.regs <- tots %>% group_by(.,pgr_damp,Centrality) %>%
    do(fit = lm((100*freq) ~ as.numeric(value), data=.)) %>%
    tidy(fit) %>% filter(.,pgr_damp == max.r)


overall.regs$term[overall.regs$term == "as.numeric(quant)"] <- "Coefficient"

overall.regs[,c(3,4,5,6)] <- sapply(overall.regs[,c(3,4,5,6)], FUN = formatC,
                                    format = "e", digits = 2)


print(xtable(overall.regs), include.rownames = F, digits = 2,
      display=c("s","s", "s", "s","g"),
      math.style.exponents = T)


tot2  <- tots %>% filter(.,Centrality %in% c("pgr.und.quant","pgr.ssc.quant"))

overall.cors <- tots %>% group_by(.,pgr_damp,Centrality) %>%
    do(fit = lm(freq ~ as.numeric(value), data=.)) %>% 
    glance(fit) 


p.difs <- data.frame(pgr_damps =  unique(overall.cors$pgr_damp), pvals = NA)
for(i in unique(overall.cors$pgr_damp)){
    r1 <- overall.cors[overall.cors$pgr_damp == i & 
                           overall.cors$Centrality =="pgr.ssc.quant",]$r.squared
    r2 <- overall.cors[overall.cors$pgr_damp == i & 
                           overall.cors$Centrality =="pgr.und.quant",]$r.squared
    out <- psych::r.test(100,(r1),(r2),twotailed = F)
    p.difs[p.difs$pgr_damps ==i,2] <- out$p
    print(out$p)
}




p2 <- ggplot(p.difs, aes(y = -log(pvals), x= p.difs[,1])) + 
    geom_line( size = 2) +
    geom_hline(yintercept=-log(0.05), linetype="dashed", color = "red") +
    labs(x = "PageRank dampening factor (alpha)", 
         y = "Negative log p-value")+
    theme_bw()+
    theme(strip.text = element_text(face="bold", size=20),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.title=element_text(face = "bold", size = 20),
          axis.text.y=element_text(size = 15),
          axis.text.x=element_text(size = 15),
          axis.ticks.y=element_blank())+ 
    annotation_custom(grobB)


library(gridExtra)
p3 <- arrangeGrob(p1,p2, nrow =1)
ggsave("images/Human_pgr_Sensitivity.pdf",p3,
       width = 18, height = 10, units = c("in"))



difs <- tot2 %>% group_by(.,pgr_damp) %>%
    do(fit = lm((100*freq) ~ as.numeric(value)/Centrality, data=.)) %>% 
    tidy(fit) 


# z <- lm((100*freq) ~ as.numeric(value), data=tot2)
# 
# 
# m.interaction <- lm((100*freq) ~ as.numeric(value)/Centrality, data=tot2)
# summary(m.interaction)
### Functionalize the above code at some point. Use the chunks below















overall.cors     <- overall.cors[,c(2,3,4,7)]
adjusted.pval    <- p.adjust(unlist(overall.cors[,4]),method = "fdr")
adjusted.pval    <- formatC(adjusted.pval, format = "e", digits = 2)
overall.cors[,4] <- formatC(unlist(pull(overall.cors[,4])), format = "e",
                            digits = 2)


print(xtable(overall.cors), include.rownames = F)

text.vals <- paste("Adj r-squared:",formatC(pull(overall.cors[,2]),digits = 2),
                   "\n p-value:", pull(overall.cors[,4]), "\n Adj p-value:",
                   adjusted.pval)

overall.cors$label <- text.vals


total <- tots %>% filter(.,pgr_damp == 0.58)




