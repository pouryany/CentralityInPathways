### There is something wrong with this test that I cannot tell.
rm(list = ls())
require(dplyr)
require(ggplot2)
require(Hmisc)
require(tidyr)
require(broom)
require(xtable)
gene.essential <- readRDS("mouse_data/MMgene_essentials.rds")




#### Getting quantiles across the data

feature.list <-  grep("norm",names(gene.essential), value = T)
feature.list <- feature.list[-c(4,6,7,8,9)]



feature.list <- gsub(".norm","",feature.list)
feature.list <- gsub("[.]","-",feature.list)
feature.list <- gsub("source","Source",feature.list)
feature.list <- gsub("sink","Sink",feature.list)
feature.list <- gsub("(^[[:alpha:]])", "\\U\\1", feature.list, perl=TRUE)
feature.list <- gsub("ssc","SSC",feature.list)
feature.list <- gsub("Pgr","PageRank",feature.list)


temp.zz <- names(gene.essential)


temp.zz <- gsub(".norm","",temp.zz)
temp.zz <- gsub("[.]","-",temp.zz)
temp.zz <- gsub("source","Source",temp.zz)
temp.zz <- gsub("sink","Sink",temp.zz)
temp.zz <- gsub("(^[[:alpha:]])", "\\U\\1", temp.zz, perl=TRUE)
temp.zz <- gsub("ssc","SSC",temp.zz)
temp.zz <- gsub("Pgr","PageRank",temp.zz)

names(gene.essential) <- temp.zz


aa <-     gene.essential %>%
    pivot_longer(., feature.list, names_to  = "Centrality", 
                 values_to  = "cent_value") %>% 
    group_by(Centrality) %>%
     mutate(.,zz = as.factor(cut2(cent_value,m = 3, g = 100))) %>%
    group_by(Centrality) %>% mutate(., zz = as.factor(zz))

    aa <- lapply(gene.essential[,feature.list],
                 FUN = function(X){as.factor(cut2(X,m=3,g=100))})
    aa <- sapply(aa, function(X){levels(X) <- 1:100; return(X)})
    aa <- cbind(aa,gene.essential$Description)
    aa <- tibble::as_data_frame(aa)
    aa2 <- tibble::as_data_frame(aa)
    aa <- tidyr::gather(aa, key = "Centrality", value = "cent_value",-c(12))

    aa$cent_value <- as.numeric(aa$cent_value)


    feat.list <- colnames(aa2)

feat.list <- feat.list[-12]
ks.list <- list()
for (i in feat.list){
    temp <- aa2[,i]
    temp <- as.numeric(pull(temp))

    this.test <- ks.test(temp[aa2$V12 == "Lethal"],
                         temp[aa2$V12 == "Normal"],
                         alternative = "less")
    ks.list <- rbind(ks.list,this.test$p.value)
}

ks.list <- data.frame(cbind(feat.list,ks.list))

ks.pvals <- formatC(unlist(ks.list$V2), digits = 3)
ks.pvals <- paste("KS p-value =",ks.pvals)


aa3 <- rbind(aa,mutate(aa,V12= "All Genes"))

ggplot(aa3,aes(cent_value,color = V12)) +
    stat_ecdf(geom = "line", size = 1.3,alpha=0.9) + facet_wrap(~Centrality ,ncol = 3)+theme_bw()+
    labs(x = "Quantile of the Z-normliazed Centrality Scores (Eq. 23)", y = "Cumulative Density")+
    guides(colour = guide_legend(override.aes = list(size=15),show.legend = "line")) +
    theme(strip.text = element_text(face="bold", size=20),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 35),
          legend.text = element_text(size = 25),
          legend.title=element_blank(),
          axis.text.y=element_text(size = 20),
          axis.text.x=element_text(size = 20),
          axis.ticks.y=element_blank(),
          legend.position=c(0.85, 0.1), legend.box = "horizontal") +
    geom_text(data=data.frame(x=50, y=0.9, label=ks.pvals,
                              Centrality = unlist(ks.list$feat.list)),
              aes(x,y,label=label),size = 7, inherit.aes=FALSE)


    ggsave("images/Mouse_KSstats.pdf",
           width = 18, height = 20, units = c("in"))




colnames(aa)



## You want a measure to say generally something works better on cancer genes.
## If you had a measure that could measure the area under CDF curve would have
## been nice.

## Remember to fix the names of pgr.source and pgr.sink

if(FALSE){
aa.cancer <-  filter(aa,  V12=="Cancer" )
aa.cancer$Centrality <- gsub(".norm","",aa.cancer$Centrality)

feat.list <- unique(aa.cancer$Centrality)
feat.list <- feat.list[-7]
ks.cancer.list <- list()
for (i in feat.list){
    temp <- aa.cancer
    temp <- as.numeric(pull(temp))

    this.test <- ks.test(temp[aa.cancer$Centrality == "pgr.ssc"],
                         temp[aa.cancer$Centrality == i],
                         alternative = "greater")
    ks.cancer.list <- rbind(ks.cancer.list,this.test$p.value)
}

ks.cancer.list  <- data.frame(cbind(feat.list,ks.cancer.list))

ks.pvals <- formatC(unlist(ks.cancer.list$V2), digits = 2)
ks.pgr.table <- cbind(feat.list,ks.pvals,ks.pvals1)
print(xtable(ks.pgr.table,digits=0),include.rownames = F)



ks.pvals <- paste("KS p-value =",ks.pvals)



?brewer.pal

library(RColorBrewer)
set.seed(2)
myColors <- brewer.pal(12,"Paired")[sample(1:12,12)]
names(myColors) <- levels(as.factor(aa.cancer$Centrality))
colScale <- scale_colour_manual(name = "grp",values = myColors)


ggplot(aa.cancer,aes(cent_value,color = Centrality)) +
    stat_ecdf(geom = "line") +
    colScale +
    labs(x = "Quantile", y = "Cumulative Density")+
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_bw()+
    theme(strip.text = element_text(face="bold", size=20),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.title=element_blank(),
          axis.text.y=element_text(size = 20),
          axis.text.x=element_text(size = 20),
          axis.ticks.y=element_blank(),
          legend.position="bottom", legend.box = "horizontal")
ggsave("images/KSstats_PGRvsAll.pdf",
       width = 18, height = 10, units = c("in"))





    aa.katz <-  filter(aa, grepl("katz",Centrality), V1=="Cancer" )

    ggplot(aa.katz,aes(cent_value,color = Centrality)) +
        stat_ecdf(geom = "point") +
        labs(x = "Quantile", y = "Cumulative Density")+
        guides(colour = guide_legend(override.aes = list(size=10))) +
        theme(strip.text = element_text(face="bold", size=20),
              plot.title = element_text(size = 20),
              axis.title = element_text(size = 25),
              legend.text = element_text(size = 20),
              legend.title=element_blank(),
              axis.text.y=element_text(size = 20),
              axis.text.x=element_text(size = 20),
              axis.ticks.y=element_blank(),
              legend.position="bottom", legend.box = "horizontal")


    ggsave("images/KSstats_katz.pdf",
           width = 18, height = 10, units = c("in"))



## Code below to be revisied.


    aa.pgr <-  filter(aa, grepl("pgr",Centrality), V1=="Cancer" )

    ggplot(aa.pgr,aes(cent_value,color = Centrality)) +
        stat_ecdf(geom = "point") +
        labs(x = "Quantile", y = "Cumulative Density")+
        guides(colour = guide_legend(override.aes = list(size=10))) +
        theme(strip.text = element_text(face="bold", size=20),
              plot.title = element_text(size = 20),
              axis.title = element_text(size = 25),
              legend.text = element_text(size = 20),
              legend.title=element_blank(),
              axis.text.y=element_text(size = 20),
              axis.text.x=element_text(size = 20),
              axis.ticks.y=element_blank(),
              legend.position="bottom", legend.box = "horizontal")
    ggsave("images/KSstats_pgr.pdf",
           width = 18, height = 10, units = c("in"))



    aa.lap <-  filter(aa, grepl("lap",Centrality), V1=="Cancer" )

    ggplot(aa.lap,aes(cent_value,color = Centrality)) +
        stat_ecdf(geom = "point") +
        labs(x = "Quantile", y = "Cumulative Density")+
        guides(colour = guide_legend(override.aes = list(size=10))) +
        theme(strip.text = element_text(face="bold", size=20),
              plot.title = element_text(size = 20),
              axis.title = element_text(size = 25),
              legend.text = element_text(size = 20),
              legend.title=element_blank(),
              axis.text.y=element_text(size = 20),
              axis.text.x=element_text(size = 20),
              axis.ticks.y=element_blank(),
              legend.position="bottom", legend.box = "horizontal")
    ggsave("images/KSstats_lap.pdf",
           width = 18, height = 10, units = c("in"))


    ggplot(gene.essential4,aes(pgr.quant2, color = Description)) +
        stat_ecdf(geom = "point", size =2) + theme_bw() + labs(x = "Quantile", y = "Cumulative Density")+
        guides(colour = guide_legend(override.aes = list(size=10))) +
        theme(strip.text = element_text(face="bold", size=20),
              plot.title = element_text(size = 20),
              axis.title = element_text(size = 25),
              legend.text = element_text(size = 20),
              legend.title=element_blank(),
              axis.text.y=element_text(size = 20),
              axis.text.x=element_text(size = 20),
              axis.ticks.y=element_blank(),
              legend.position="bottom", legend.box = "horizontal") +
        annotate("text",x=50, y= 1, size = 10,  label= "KS p-value < 4e-04")



    p.vals.process <- gene.essential
    feature.list <-  grep("norm",names(gene.essential), value = T)
    feature.list <- feature.list[-c(4,6,7,8,9)]

    p.vals.process <- p.vals.process %>%
        gather_(., key = "Centrality", value = "cent_value",feature.list)%>%
        group_by(pathway.name,Centrality) %>%
        do(pval =tidy(wilcox.test( cent_value~ Description,
                              alternative = "greater",paired = F,exact=FALSE, data = .))) %>%
        unnest()  %>%
        group_by(Centrality) %>%
        mutate(., fdr = p.adjust(p.value)) %>%
        filter(., fdr < 0.05)


    # Creating confusion matrix
    mat1 <- table(p.vals.process$Centrality,p.vals.process$pathway.name)
    mat2 <- table(p.vals.process$pathway.name,p.vals.process$Centrality)


    confusion <- as.matrix(mat1) %*% as.matrix(mat2)
    colnames(confusion) <- gsub(".norm","",colnames(confusion))
    rownames(confusion) <- colnames(confusion)

    print(xtable(confusion,digits=0,rotate.colnames = TRUE))











mutate(., fdr = p.adjust(p.value)) %>%
    filter(., fdr < 0.05)

?cut2

zz <- cut2(gene.essential$pgr.ssc.norm,m = 3, g = 100)
levels(zz) <- 1:100
zz2 <- cut2(gene.essential$katz.ssc.norm,m = 3, g = 100)
levels(zz2) <- 1:100

zz3 <- cut2(gene.essential$pgr.und.norm,m = 3, g = 100)
levels(zz3) <- 1:100


# Creating quantile scores
gene.essential <-  mutate(gene.essential, pgr.quant2 = as.numeric(zz))
gene.essential <-  mutate(gene.essential, katz.ssc.quant2 = as.numeric(zz2))
gene.essential <-  mutate(gene.essential, pgr.quant3 = as.numeric(zz3))



# Comparing the quantile scores
zzz0 <-  gene.essential %>% dplyr::select(.,katz.ssc.quant2) %>%
    unlist %>% unname #%>% ecdf()

zzz1 <-  gene.essential %>% filter(., Description =="Normal") %>%
    dplyr::select(.,katz.ssc.quant2) %>% unlist %>% unname #%>% ecdf()

zzz2<-   gene.essential %>% filter(., Description =="Cancer") %>%
    dplyr::select(.,katz.ssc.quant2) %>% unlist %>% unname #%>% ecdf()

zzz3 <-  gene.essential %>% filter(., Description =="Normal") %>%
    dplyr::select(.,pgr.quant3) %>% unlist  #%>% ecdf()

zzz4 <-  gene.essential %>% filter(., Description =="Cancer") %>%
    dplyr::select(.,pgr.quant3) %>% unlist  #%>% ecdf()


# Just a check between KS Statistics
aaaa <-ks.test(zzz2,zzz1, alternative = "less")
aaaa$statistic


# Getting KS Statistics across models
gene.essential2 <- rbind(gene.essential, mutate(gene.essential, Description = "All genes"))

ssk <- ks.test(gene.essential2$katz.ssc.quant2[gene.essential2$Description == "Cancer"],
               gene.essential2$katz.ssc.quant2[gene.essential2$Description == "Normal"],
               alternative = "less")
ssp <- ks.test(gene.essential2$pgr.quant2[gene.essential2$Description == "Cancer"],
               gene.essential2$pgr.quant2[gene.essential2$Description == "Normal"],
               alternative = "less")
udp <- ks.test(gene.essential2$pgr.quant3[gene.essential2$Description == "Cancer"],
               gene.essential2$pgr.quant3[gene.essential2$Description == "Normal"],
               alternative = "less")

ks.pvals <- c(ssp$p.value,udp$p.value,ssk$p.value)
ks.pvals <- formatC(ks.pvals, digits = 3)
ks.pvals <- paste("KS p-value <",ks.pvals)

# Getting plots out!
gene.essential2 <- gene.essential2 %>% gather(.,key = "Type", value = "ranks", pgr.quant2,pgr.quant3,katz.ssc.quant2)
gene.essential2$Type <- as.factor(gene.essential2$Type)
levels(gene.essential2$Type) <- c("Source-Sink PageRank", "Undirected PageRank", "Source-Sink Katz")

ggplot(gene.essential2,aes(ranks,color = Description)) +
    stat_ecdf(geom = "point") + facet_wrap(~Type ,ncol = 3)+theme_bw() + labs(x = "Quantile", y = "Cumulative Density")+
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme(strip.text = element_text(face="bold", size=20),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.title=element_blank(),
          axis.text.y=element_text(size = 20),
          axis.text.x=element_text(size = 20),
          axis.ticks.y=element_blank(),
          legend.position="bottom", legend.box = "horizontal") +
    geom_text(data=data.frame(x=50, y=0, label=ks.pvals,
                              Type = c("Source-Sink PageRank", "Undirected PageRank", "Source-Sink Katz")),
              aes(x,y,label=label),size = 7, inherit.aes=FALSE)



ggplot(gene.essential,aes(pgr.und.quant,pgr.source.quant)) + geom_point()



# Some more plotting

gene.essential3 <- mutate(gene.essential) %>% filter(.,Description == "Cancer") %>% mutate(., Description = "Source-Sink Cancer")
gene.essential4 <- mutate(gene.essential3, pgr.quant2 = pgr.quant3 , Description = "Undirected Cancer")
gene.essential4 <- rbind(gene.essential3, gene.essential4)

ssp.udp <- ks.test(gene.essential4$pgr.quant2[gene.essential4$Description == "Source-Sink Cancer"],
                   gene.essential4$pgr.quant2[gene.essential4$Description == "Undirected Cancer"],
                   alternative = "less")
ssp.udp$p.value

ggplot(gene.essential4,aes(pgr.quant2, color = Description)) +
    stat_ecdf(geom = "point", size =2) + theme_bw() + labs(x = "Quantile", y = "Cumulative Density")+
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme(strip.text = element_text(face="bold", size=20),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.title=element_blank(),
          axis.text.y=element_text(size = 20),
          axis.text.x=element_text(size = 20),
          axis.ticks.y=element_blank(),
          legend.position="bottom", legend.box = "horizontal") +
    annotate("text",x=50, y= 1, size = 10,  label= "KS p-value < 4e-04")

ks.test(zzz4,zzz2, alternative = "greater")
}
