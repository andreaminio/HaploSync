#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(htmlwidgets))
suppressPackageStartupMessages(library(plyr))

options(stringsAsFactors = F)

diploid_trace <- read.table("diploid_gene_count_trace.hap1.txt",sep="\t",header=T)
colnames(diploid_trace) <- c("CHR","START","STOP","ID","HAP1","HAP2","DESC","DESC_COUNT")
save.image( file = "diploid_gene_count_trace.hap1.rda")
llply(unique(diploid_trace$CHR), function(x){
  df.x <- subset(diploid_trace, CHR==x)
  p <- ggplot(data=df.x) +
    geom_step(aes(x=START,y=HAP1,color="HAP1")) +
    geom_step(aes(x=START,y=HAP2,color="HAP2")) +
    labs(title=x, x ="Gene start (bp)", y = "Hit count (iden > 95% && cov > 95%)")
  saveWidget(ggplotly(p, dynamicTicks = TRUE) , paste0(x,".html"))
})
