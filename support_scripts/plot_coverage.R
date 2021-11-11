#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(htmlwidgets))

repeats_file <- Sys.glob("*.repeat.bed.gz")
gap_file <- Sys.glob("*.gap.bed")
raw_cov_file <- Sys.glob("*.cov.bed.gz")
smoothed_cov_file <- Sys.glob("*.smooth_cov.bed.gz")

repeats_bed <- read.table( gzfile(repeats_file), header=F)
colnames(repeats_bed) <- c("Chr","Start","Stop")
gap_bed <- read.table( gap_file , header=F)
colnames(gap_bed) <- c("Chr","Start","Stop")
raw_cov_trace <- read.table( gzfile(raw_cov_file) , header=F)
colnames(raw_cov_trace)<- c("Chr","Start","Stop","Coverage" )
smoothed_cov_trace <- read.table( gzfile(smoothed_cov_file) , header=F)
colnames( smoothed_cov_trace ) <- c("Chr","Start","Stop","Coverage" )

max_y = max( max(raw_cov_trace$Coverage) , max(smoothed_cov_trace$Coverage) ) + 100

p <- ggplot() +
  geom_rect(data=repeats_bed, aes(ymin=-10, ymax=max_y, xmin=Start, xmax=Stop , fill="Repeats" ) , alpha =0.5) +
  geom_rect(data=gap_bed    , aes(ymin=-10, ymax=max_y, xmin=Start, xmax=Stop , fill="Gaps" ) , alpha =0.5) +
  geom_step(data=raw_cov_trace      , aes (x=Start , y=Coverage , colour = "Raw coverage") ) +
  geom_step(data=smoothed_cov_trace , aes (x=Start , y=Coverage , colour = "Smoothed coverage") ) +
  labs(title="Sequencig data coverage", x ="Genome position (bp)", y = "Coverage")

save.image( file = "coverage_plot.rda" )

saveWidget(ggplotly(p, dynamicTicks = TRUE) , "coverage_plot.html")
