### Plotting correlations and distributions

require(dplyr)
require(ggplot2)

gene.essential <- readRDS("Human Processing/gene_essentials.rds")


norm.names <-  names(gene.essential)[grep("norm",names(gene.essential))]

for (i in norm.names){
    this.plot <- ggplot(gene.essential,aes_string(c(i),   colour = "pathway.name")) +
        stat_density( kernel = "g")+
        theme_bw() +
        theme(strip.text = element_text(face="bold", size=20),
              plot.title = element_text(size = 20),
              axis.title = element_text(size = 25),
              legend.position="none",
              axis.text.y=element_text(size = 20),
              axis.text.x=element_text(size = 20),
              axis.ticks.y=element_blank())
    file.name <- paste0("images/normalized/",gsub("\\.","_", i), ".pdf")
    ggsave(file.name,this.plot, width = 18, height = 10, units = c("in"))
}



vec.names <-  names(gene.essential)[grep("vec",names(gene.essential))]

for (i in vec.names){
    this.plot <- ggplot(gene.essential,aes_string(c(i),   colour = "pathway.name")) +
        stat_density( kernel = "g")+
        theme_bw() +
        theme(strip.text = element_text(face="bold", size=20),
              plot.title = element_text(size = 20),
              axis.title = element_text(size = 25),
              legend.position="none",
              axis.text.y=element_text(size = 20),
              axis.text.x=element_text(size = 20),
              axis.ticks.y=element_blank())
    file.name <- paste0("images/vectors/",gsub("\\.","_", i), ".pdf")
    ggsave(file.name,this.plot, width = 18, height = 10, units = c("in"))
}


quant.names <-  names(gene.essential)[grep("quant",names(gene.essential))]

for (i in quant.names){
    this.plot <- ggplot(gene.essential,aes_string(c(i),   colour = "pathway.name")) +
        stat_density( kernel = "g")+
        theme_bw() +
        theme(strip.text = element_text(face="bold", size=20),
              plot.title = element_text(size = 20),
              axis.title = element_text(size = 25),
              legend.position="none",
              axis.text.y=element_text(size = 20),
              axis.text.x=element_text(size = 20),
              axis.ticks.y=element_blank())
    file.name <- paste0("images/quants/",gsub("\\.","_", i), ".pdf")
    ggsave(file.name,this.plot, width = 18, height = 10, units = c("in"))
}



# Plotting Correlations
#biocLite("GGally")
library(GGally)
gene.cor.data <- gene.essential %>% 
                 dplyr::select(., deg.quant,
                               katz.sink.quant,katz.source.quant,katz.ssc.quant,
                               lap.sink.quant,lap.source.quant ,lap.ssc.quant, lap.und.quant,
                               pgr.sink.quant,pgr.source.quant, pgr.ssc.quant, pgr.und.quant)                               
ggpairs(gene.cor.data, method = "pearson", columns = 1:ncol(gene.cor.data),
        columnLabels = c("Degree", "Katz-Sink", "Katz-Source","Katz-SSC",
                         "Lap-Sink","Lap-Source","Lap-SSC","Lap-Und",
                         "Pgr-Source","Pgr-Sink", "Pgr-SSC", "Pgr-Und"),
        title = "",
        axisLabels = "internal",
        diag  = list(continuous=wrap("text", labelSize=30)),
        upper = list(continuous=wrap(ggally_cor,size=5)),
        lower = list(continuous = wrap("points", alpha = 0.3,    size=0.1))
)  + theme_bw()
ggsave("images/corrsALL.pdf",
       width = 18, height = 10, units = c("in"))






