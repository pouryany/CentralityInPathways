rm(list = ls())
library(dplyr)
library(broom)
library(ggplot2)
library(xtable)

gene.essential <- readRDS("Mouse Processing/MMgene_essentials.rds")








gene.essential$Description <- factor(gene.essential$Description)
# levels(gene.temp$Description)
 gene.essential$Description <- relevel(gene.essential$Description,ref = "Normal")


## The regression analysis

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

### Functionalize the above code at some point. Use the chunks below


var  <- "katz.ssc.quant"
cent <- "Katz_Source_Sink"
Desc <- "Description"

quant.sorter <- function(data,var, Desc,cent){

    varval1 <- lazyeval::interp(~quant, quant = as.name(var))
    varval2 <- lazyeval::interp(~Centrality, Centrality = cent)
    sorted  <- gene.essential %>% group_by_(var, as.name(Desc))  %>%
        summarise(n = n()) %>% mutate(freq = n/ sum(n)) %>%
        filter_(., ~Description == "Lethal") %>%
        mutate_(.dots = setNames(list(varval1,varval2),
                                 c("quant", "Centrality")))

    return(sorted)
}


model.vars <-  c("katz.ssc.quant", "katz.source.quant","katz.sink.quant",
                 "deg.quant","beet.source.quant","beet.sink.quant",
                 "beet.ssc.quant","beet.und.quant","pgr.source.quant",
                 "pgr.sink.quant","pgr.ssc.quant","pgr.und.quant",
                 "lap.source.quant","lap.sink.quant","lap.ssc.quant",
                 "lap.und.quant")

model.names <-  c("katz_ssc", "katz_source","katz_sink",
                  "deg","beet_source","beet_sink",
                  "beet_ssc","beet_und","pgr_source",
                  "pgr_sink","pgr_ssc","pgr_und",
                  "lap_source","lap_sink","lap_ssc",
                  "lap_und")







feature.list <-  grep("quant",names(gene.essential))
-c(1,2,3,4,58,59,60,61)

reg.table <- gene.essential %>%
    gather(., key = "Centrality", value = "cent_value",
           gather_cols=feature.list) %>%
    group_by(Centrality) %>%
    do(fit =tidy(glm(as.factor(Description) ~ cent_value, family="binomial", data = .))) %>%
    unnest()  %>%
    group_by(Centrality)


mutate(., fdr = p.adjust(p.value)) %>%
    filter(., fdr < 0.05)


bbb <- reg.table %>% filter(., term != "(Intercept)") %>%
    mutate(., estimate = exp(estimate))
levels(gene.essential$Description) <- c(0,1)

glm(as.factor(Description) ~ katz.sink.norm + pgr.sink.norm, family= binomial, data = gene.essential)



log(log(log(katz.ssc.vec)) +2)
as.logical(factor(gene.essential$Description))
log(100)
exp(0.37)
ggplot(gene.essential,aes(x =katz.ssc.vec ,   colour = Description)) +
    geom_density( kernel = "g")+
    theme_bw()

ggplot(gene.essential,aes(x =log(katz.source.vec) , y = log(katz.sink.vec),
                          colour = Description)) + geom_point()
ggplot(gene.essential,aes(x =rgamma(nrow(gene.essential),3.27839639,1.00738148 ))) +
    geom_density( kernel = "g")


?dgamma

library(MASS)
my.mle<-fitdistr(gene.essential$katz.ssc.vec, densfun="gamma")
my.mle$sd
BIC(my.mle)


feature.list <-  grep("quant",names(gene.essential),value = T)
feature.list <- feature.list[-2]
reg.res <- gene.essential %>% mutate(., deg.quant = all.degree)
reg.res <-    gene.essential %>%
    gather_(., key = "Centrality", value = "cent_value",
            gather_cols=feature.list)


reg.res %>% group_by(Centrality) %>% group_by(Description,cent_value) %>%
    summarise(n = n()) %>%
    mutate(freq = n/ sum(n)) %>%
    filter(., Description== "Lethal") %>% mutate(.,Models = Centrality)


reg.res %>% group_by(., Centrality) %>%
    do(fit = lm(freq ~ as.numeric(cent_value), data=.)) %>%    glance(fit)


