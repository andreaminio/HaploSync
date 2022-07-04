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

min_len <- as.numeric(params$min_align)
min_iden <- as.numeric(params$similarity)

structure_file <- params$structure
legacy_file <- params$legacy
markers_all_file <- params$markers
markers_dup_file <- params$dup_markers

# load data
tb_1 <- read_delim(file = nucmer_file,
                 delim = "\t",
                 col_names = T)
coords1 <- tb_1 %>% filter(tID==nucmer_ref , qID==nucmer_query , identity > min_iden , match > min_len )

if (nucmer_file_self == "") {
  coords2 <- 0
} else {
  tb_2 <- read_delim(file = nucmer_file_self,
                 delim = "\t",
                 col_names = T)
  coords2 <- tb_2 %>% filter(tID==nucmer_ref , qID==nucmer_ref , identity > min_iden , match > min_len )
}
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

all_marker_trace <- read_delim(markers_all_file,
                            delim = "\t",
                            col_names = c( "Chr" , "Start" , "Stop" , "Marker_Id") ,
                            )
dup_markers_trace <- read_delim(markers_dup_file,
                            delim = "\t",
                            col_names = c( "Chr" , "Start" , "Stop" , "Marker_Id") ,
                            )




# subset
refLen <- max(coords1$tLen)
queryLen <- max(coords1$qLen)

all_marker_trace_r <- all_marker_trace %>% filter(Chr==nucmer_ref )
if (dim(dup_markers_trace)[1] == 0) {
  dup_markers_trace_r <- dup_markers_trace
} else{
   dup_markers_trace_r <- dup_markers_trace %>% filter(Chr==nucmer_ref )
}


# plot
p0 <- ggplot(coords1) + geom_segment(
    aes( x = tStart,
         y = qStart,
         xend = tStop,
         yend = qStop,
         color = identity ,
        )
  ) +
  scale_color_viridis_c(name = "%id", direction = -1, begin = 0) +
  scale_x_continuous(limits = c(0, refLen)) +
  scale_y_continuous(limits = c(0, max(refLen , queryLen) )) +
  labs(title = paste0( nucmer_ref , " - Self alignment") , x = "" , y = nucmer_ref)

if (!(coords2==0)) {
p0 <- p0 + geom_segment(data = coords2 ,
     aes( x = tStart,
         y = qStart,
         xend = tStop,
         yend = qStop,
         color = identity ,
        )
  )
}

p0 <- p0 + theme(plot.title = element_text(color="black", size=8, face="plain", margin=margin(t=1,r=0,b=1,l=0,"mm")),
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
gt0 <- ggplotGrob(p0)
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

gt0_widths <- gt0$widths
gt5_widths <- gt5$widths
gt6_widths <- gt6$widths

maxWidth <- unit.pmax(gt0_widths, gt5_widths , gt6_widths )

gt0$widths <- maxWidth
gt5$widths <- maxWidth
gt6$widths <- maxWidth

layout <- rbind(1,2,3)

g <- arrangeGrob( gt0, gt5 , gt6 , layout_matrix=layout , heights=c(.8, .1 , .1 ) )
ggsave(file=paste0(out_dir,"/",out_file), g , height = 12 , width = 10 )