# load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(optparse))

# prepare variables

option_list <- list(
  make_option(c("-c","--coords"), type="character", default=NULL,
              help="Mapping coords table for different genotyopes [default %default]",
              dest="coords"),
  make_option(c("-s","--self"), type="character", default=NULL,
              help="Coordinates of target sequence self alignment [default %default]",
              dest="coords_self"),
  make_option(c("-o","--out"), type="character", default="out",
              help="Output file name [default %default]",
              dest="out_file"),
  make_option(c("-d","--output_dir"), type="character", default=".",
              help="Output directory [default %default]",
              dest="out_dir"),
  make_option(c("-t", "--target"), type="character", default=NULL,
              help="Reference sequence ID [default %default]",
              dest="refID"),
  make_option(c("-q", "--query"), type="character", default=NULL,
              help="Query sequence ID [default %default]",
              dest="queryID"),
  make_option(c("-m", "--min-alignment-length"), type="numeric", default=3000,
              help="Minimum hit length cutoff in bp [default %default]",
              dest="min_align"),
  make_option(c("-i","--identity"), type="numeric", default=90,
              help="minimum idetity of hits to plot in percentage [default %default]",
              dest="similarity"),
  make_option(c("-g","--genecounts"), type="character", default=NULL,
              help="Gene mapping counts table [default %default]",
              dest="counts"),
  make_option(c("-r", "--ratio"), type="numeric", default=0.33,
              help="Ratio between haplotypes to call a gene impared (range 0-1) [default %default]",
              dest="ratio"),
  make_option(c("-a", "--structure"), type="character", default=NULL,
              help="Target sequence structure [default %default]",
              dest="structure"),
  make_option(c("-l", "--legacy"), type="character", default=NULL,
              help="Target legacy sequence structure [default %default]",
              dest="legacy"),
  make_option(c("-b", "--markers"), type="character", default=NULL,
              help="Position of markers on target sequence [default %default]",
              dest="markers"),
  make_option(c("-e", "--duplicated_markers"), type="character", default=NULL,
              help="Position of duplicated markers on target sequence [default %default]",
              dest="dup_markers")
)


options(error=traceback)
parser <- OptionParser(usage = "%prog -c alignments.coords.txt -s selfmap.coords.txt -o out [options]",option_list=option_list)
params <- parse_args(parser)

out_file <- params$out_file
out_dir <- params$out_dir
nucmer_file <- params$coords
nucmer_file_self <- params$coords_self
nucmer_query <- paste0(params$queryID)
nucmer_ref <- paste0(params$refID)
nucmer_title <- paste0(params$queryID," on " , params$refID)

counts_file <- params$counts
gmap_hap1 <- params$refID

min_len <- as.numeric(params$min_align)
min_iden <- as.numeric(params$similarity)
ratio <- as.numeric(params$ratio)

structure_file <- params$structure
legacy_file <- params$legacy
markers_all_file <- params$markers
markers_dup_file <- params$dup_markers

category_conversion_1= c(
"dup_1_0" = "CDS with no hit in Hap2" ,
"dup_gt1_0" = "CDS with multiple hits in Hap1 but no hit in Hap2" ,
"dup_gt1_gt0_grt" = paste0("CDS present in both haplotypes, Hap1/Hap2 counts ratio > ",(1/ratio)) ,
"dup_gt1_gt0_lrt" = paste0("CDS present in both haplotypes, Hap1/Hap2 counts ratio < ",ratio)
)

category_conversion_2= c(
"dup_0_1" = "CDS with no hit in Hap1" ,
"dup_0_gt1" = "CDS with multiple hits in Hap2 but no hit in Hap1" ,
"dup_gt0_gt1_grt" = paste0("CDS present in both haplotypes, Hap2/Hap1 counts ratio > ",(1/ratio)) ,
"dup_gt0_gt1_lrt" = paste0("CDS present in both haplotypes, Hap2/Hap1 counts ratio < ",ratio)
)

# load data
tb <- read_delim(file = nucmer_file,
                 delim = "\t",
                 col_names = T)
# Columns:
# tID
# tLen
# tStart
# tStop
# qID
# qLen
# qStart
# qStop
# identity
# match

diploid_trace_hap1 <- read_delim(counts_file,
                            delim = "\t",
                            col_names = T)

# Columns:
# Chr
# Start
# Stop
# Gene_Id
# Hap1_count
# Hap2_count
# Hap2_to_Hap1_ratio
# Description
# Gene_count_with_description

all_marker_trace <- read_delim(markers_all_file,
                            delim = "\t",
                            col_names = c( "Chr" , "Start" , "Stop" , "Marker_Id") ,
                            )
dup_markers_trace <- read_delim(markers_dup_file,
                            delim = "\t",
                            col_names = c( "Chr" , "Start" , "Stop" , "Marker_Id") ,
                            )



# subset
coords2 <- tb %>% filter(tID==nucmer_ref , qID==nucmer_query , identity > min_iden , match > min_len )

refLen <- max(coords2$tLen)
queryLen <- max(coords2$qLen)

diploid_trace_hap1 <- diploid_trace_hap1 %>%
  mutate(test = Hap1_count/Hap2_count) %>%
  mutate(dup_1_0 = (Hap1_count==1 & Hap2_count==0),
         dup_gt1_0 = (Hap1_count>1 & Hap2_count==0),
         dup_gt1_gt0_grt = (Hap1_count>1 & Hap2_count>0 & test > (1/ratio) ) ,
         dup_gt1_gt0_lrt = (Hap1_count>1 & Hap2_count>0 & test < ratio ) ) %>%
  filter(Chr==gmap_hap1)

diploid_trace2 <- pivot_longer(diploid_trace_hap1, cols = c("Hap1_count", "Hap2_count"), names_to = "HAP", values_to = "count")

diploid_trace3 <- pivot_longer(diploid_trace_hap1, cols = c("dup_1_0", "dup_gt1_0", "dup_gt1_gt0_grt" , "dup_gt1_gt0_lrt" ), names_to = "DUP", values_to = "DUP_test")

diploid_trace3 <- diploid_trace3 %>%
  mutate(DUP_test = str_replace_all(DUP_test, "FALSE", "0")) %>%
  mutate(DUP_test = str_replace_all(DUP_test, "TRUE", "1")) %>%
  mutate(DUP_test = as.numeric(DUP_test))

all_marker_trace_r <- all_marker_trace %>% filter(Chr==nucmer_ref )
dup_markers_trace_r <- dup_markers_trace %>% filter(Chr==nucmer_ref )


# plot


p1 <- ggplot(coords2) + geom_segment(
    aes( x = tStart,
         y = qStart,
         xend = tStop,
         yend = qStop,
         color = identity ,
        )
  ) +
  scale_color_viridis_c(name = "%id", direction = -1, begin = 0) +
  scale_x_continuous(limits = c(0, refLen)) +
  scale_y_continuous(limits = c(0, queryLen)) +
  labs(title = nucmer_title, x = paste0( nucmer_ref , " - Sequence position (bp)")  , y = nucmer_query)

p1 <- p1 + theme(plot.title = element_text(color="black", size=8, face="plain", margin=margin(t=1,r=0,b=1,l=0,"mm")),
               axis.title.x = element_text(color="black", size=7, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm")),
               axis.title.y = element_text(color="black", size=7, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm")),
               axis.text.x = element_text(color="black", size=8, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm"), hjust = 0.5),
               axis.text.y = element_text(color="black", size=8, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm"), vjust = 0.5),
               axis.line.x = element_line(color="black", size=0.5, linetype="solid"),
               axis.line.y = element_line(color="black", size=0.5, linetype="solid"),
               axis.ticks.x = element_line(color="black", size=0.5, linetype="solid"),
               axis.ticks.y = element_line(color="black", size=0.5, linetype="solid"),
               plot.background = element_blank(),
               panel.background = element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.grid.minor.x = element_blank(),
               plot.margin=unit(c(1,1,1,1),"mm"),
               legend.position="right")


p2 <- ggplot(diploid_trace2,
    aes(
      x = Start,
      y = count,
      color = HAP,
      group = 1,
    )
  ) +
  geom_step() +
  scale_colour_manual(values = brewer.pal(8, "Dark2")[1:2] , labels = c("Hap 1 counts" , "Hap 1 counts") ) +
  scale_y_continuous(breaks = seq(min(diploid_trace2$count), max(diploid_trace2$count), 1)) +
  labs(title = paste0( nucmer_ref , " - Gene mapping counts")  , x = paste0( nucmer_ref , " - Sequence position (bp)"), y = "Gene alignments counts" , color = "")

p2 <- p2 +
          theme(plot.title = element_text(color="black", size=8, face="plain", margin=margin(t=1,r=0,b=1,l=0,"mm")),
                 axis.title.x = element_text(color="black", size=7, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm")),
                 axis.title.y = element_text(color="black", size=7, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm")),
                 axis.text.x = element_text(color="black", size=8, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm"), hjust = 0.5),
                 axis.text.y = element_text(color="black", size=8, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm"), vjust = 0.5),
                 axis.line.x = element_line(color="black", size=0.5, linetype="solid"),
                 axis.line.y = element_line(color="black", size=0.5, linetype="solid"),
                 axis.ticks.x = element_line(color="black", size=0.5, linetype="solid"),
                 axis.ticks.y = element_line(color="black", size=0.5, linetype="solid"),
                 plot.background = element_blank(),
                 panel.background = element_blank(),
                 panel.grid.major.y = element_blank(),
                 panel.grid.major.x = element_blank(),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 plot.margin=unit(c(2,2,2,2),"mm"),
                 legend.position="right"
            )


#- dup_1_0: CDS of hap1 absent in hap2
#- dup_gt1_0: CDS present more than once in hap1 but absent in hap2
#- dup_gt1_gt0_rtgt2: CDS present more than once in hap1, at least once in hap2, with a ratio hap1/hap2 greater than 2


p3 <- ggplot(diploid_trace3 ,
    aes(
    x = Start,
    y = DUP_test,
    color = DUP,
    group = 1,
    )
  ) +
  geom_step() +
  facet_wrap(~DUP , ncol=1 , scales = "free_y", labeller = as_labeller(category_conversion_1) ) +
  scale_y_continuous(breaks = c(0,1) , limits=c(0,1.1) , labels = c("False" , "True")) +
  scale_colour_manual(values = brewer.pal(8, "Dark2")[3:6] ) +
  labs(x = paste0( nucmer_ref , " - Sequence position (bp)"), y = "Test result" , color = "")

p3 <- p3 +
        theme(plot.title = element_text(color="black", size=8, face="plain", margin=margin(t=1,r=0,b=1,l=0,"mm")),
                 axis.title.x = element_text(color="black", size=7, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm")),
                 axis.title.y = element_text(color="black", size=7, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm")),
                 axis.text.x = element_text(color="black", size=8, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm"), hjust = 0.5),
                 axis.text.y = element_text(color="black", size=8, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm"), vjust = 0.5),
                 axis.line.x = element_line(color="black", size=0.5, linetype="solid"),
                 axis.line.y = element_line(color="black", size=0.5, linetype="solid"),
                 axis.ticks.x = element_line(color="black", size=0.5, linetype="solid"),
                 axis.ticks.y = element_line(color="black", size=0.5, linetype="solid"),
                 plot.background = element_blank(),
                 panel.background = element_rect(fill = "white", size = 0.5, linetype = "solid"),
                 panel.grid.major.y = element_line(color="#777777", size=0.2, linetype="solid"),
                 panel.grid.major.x = element_blank(),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 plot.margin=unit(c(2,2,2,2),"mm"),
                 legend.position="none",
                 panel.spacing = unit(.8, "lines") ,
                 strip.background = element_rect(fill = "#999999" ),
                 strip.text.x = element_text(margin = margin(.3, 0, .1, 0, "lines") , color = "white ")
    )


p5 <- ggplot() +
  geom_rect( data=all_marker_trace_r ,
    aes(
      xmin = Start,
      xmax = Stop,
      ymin = 0 ,
      ymax = 1 ,
      colour = "Marker",
      fill = "Marker"
        ) ,
    show.legend = FALSE
  ) +
  scale_color_manual(name = "", breaks = c( "Marker" , "Duplicated marker" ) , values = c("black" ,"red") ) +
  scale_fill_manual(name = "", breaks = c( "Marker" , "Duplicated marker" ) , values = c("black" ,"red") ) +
  scale_x_continuous(limits = c(0, refLen)) +
  scale_y_continuous(breaks = c( 0 , 1)) +
  labs(title = "Markers" , x = "Sequence position (bp)", y = "")

p5 <- p5 +
          theme(plot.title = element_text(color="black", size=8, face="plain", margin=margin(t=1,r=0,b=1,l=0,"mm")),
                 axis.title.x = element_text(color="black", size=8, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm")),
                 axis.title.y = element_text(color="black", size=8, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm")),
                 axis.text.x = element_text(color="black", size=6, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm"), hjust = 0.5),
                 axis.text.y = element_text(color="black", size=6, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm"), vjust = 0.5),
                 axis.line.x = element_line(color="black", size=0.5, linetype="solid"),
                 axis.line.y = element_line(color="black", size=0.5, linetype="solid"),
                 axis.ticks.x = element_line(color="black", size=0.5, linetype="solid"),
                 axis.ticks.y = element_line(color="black", size=0.5, linetype="solid"),
                 plot.background = element_blank(),
                 panel.background = element_blank(),
                 panel.grid.major.y = element_line(color="gray80", size=0.2, linetype="solid") ,
                 panel.grid.major.x = element_line(color="gray80", size=0.2, linetype="solid"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 plot.margin=unit(c(2,2,2,2),"mm"),
                 #legend.position="right"
            )

if (dim(dup_markers_trace_r)[1] == 0) {
  filler = tibble(Start = 1 , Stop = refLen)
  p6 <- ggplot() +
    geom_rect( data= filler,
               aes(
                 xmin = Start,
                 xmax = Stop,
                 ymin = 0 ,
                 ymax = 1 ,
                 fill = "No marker",
               ) ,
              show.legend = FALSE
    )

} else {
  p6 <- ggplot() +
    geom_rect( data=dup_markers_trace_r ,
               aes(
                 xmin = Start,
                 xmax = Stop,
                 ymin = 0 ,
                 ymax = 1 ,
                 colour = "Duplicated marker",
                 fill = "Duplicated marker",
                ),
               show.legend = FALSE
    )
}

p6 <- p6 +
  scale_color_manual(name = "", breaks = c( "Duplicated marker" ) , values = c("red") ) +
  scale_fill_manual(name = "", breaks = c( "Duplicated marker" , "No marker" ) , values = c("red" , "gray95") ) +
  scale_x_continuous(limits = c(0, refLen)) +
  scale_y_continuous(breaks = c( 0 , 1)) +
  labs(title = "Duplictaed markers" , x = "", y = "")

p6 <- p6 +
          theme(plot.title = element_text(color="black", size=8, face="plain", margin=margin(t=1,r=0,b=1,l=0,"mm")),
                 axis.title.x = element_text(color="black", size=8, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm")),
                 axis.title.y = element_text(color="black", size=8, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm")),
                 axis.text.x = element_text(color="black", size=6, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm"), hjust = 0.5),
                 axis.text.y = element_text(color="black", size=6, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm"), vjust = 0.5),
                 axis.line.x = element_line(color="black", size=0.5, linetype="solid"),
                 axis.line.y = element_line(color="black", size=0.5, linetype="solid"),
                 axis.ticks.x = element_line(color="black", size=0.5, linetype="solid"),
                 axis.ticks.y = element_line(color="black", size=0.5, linetype="solid"),
                 plot.background = element_blank(),
                 panel.background = element_blank(),
                 panel.grid.major.y = element_line(color="gray80", size=0.2, linetype="solid") ,
                 panel.grid.major.x = element_line(color="gray80", size=0.2, linetype="solid"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 plot.margin=unit(c(2,2,2,2),"mm"),
                 #legend.position="right"
            )


pdf(NULL)
gt1 <- ggplotGrob(p1)
gt2 <- ggplotGrob(p2)
gt3 <- ggplotGrob(p3)
gt5 <- ggplotGrob(p5)
gt6 <- ggplotGrob(p6)

#pdf("/Users/andreaminio/Desktop/tmp_transfer/test_ChrBoard.pdf", height = 20 , width = 10 )
#newWidth = unit.pmax(gt1$widths[2:3], gt2$widths[2:3] ,gt2$widths[2:3], gt3$widths[2:3] )
#gt0$widths[2:3] = as.list(newWidth)
#gt1$widths[2:3] = as.list(newWidth)
#gt2$widths[2:3] = as.list(newWidth)
#gt3$widths[2:3] = as.list(newWidth)
#grid.draw(rbind(gt0, gt1 , gt2 , gt3))
#dev.off()

gt1_widths <- gt1$widths
gt2_widths <- gt2$widths
gt3_widths <- gt3$widths
gt5_widths <- gt5$widths
gt6_widths <- gt6$widths

maxWidth <- unit.pmax( gt1_widths, gt2_widths, gt3_widths , gt5_widths , gt6_widths )

gt1$widths <- maxWidth
gt2$widths <- maxWidth
gt3$widths <- maxWidth
gt5$widths <- maxWidth
gt6$widths <- maxWidth

layout <- rbind(1,2,3,4,5,6,7,8)

g <- arrangeGrob( gt1 , gt5 , gt6 , gt2 , gt3, layout_matrix=layout , heights=c(.25, .05 , .05 , .1 , .2) )
ggsave(file=paste0(out_dir,"/",out_file), g , height = 20 , width = 10 )