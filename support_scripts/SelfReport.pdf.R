# load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plyr))

# prepare variables
# " -d " + output_dir + " -o " + output_file + "-c " + coords + " -g " + gene_counts + " -m " + markers_all + " -n " + markers_dup + " -q " + queryID + " -s " + structure
option_list <- list(
  make_option(c("-o","--out"), type="character", default="out",
              help="Output file name [default %default]",
              dest="out_file"),
  make_option(c("-d","--output_dir"), type="character", default=".",
              help="Output directory [default %default]",
              dest="out_dir"),
  make_option(c("-c","--coords"), type="character", default=NULL,
              help="Mapping coords table for different genotyopes [default %default]",
              dest="coords"),
  make_option(c("-g","--genecounts"), type="character", default=NULL,
              help="Gene mapping counts table [default %default]",
              dest="gene_counts"),
  make_option(c("-m", "--markers"), type="character", default=NULL,
              help="Markers position [default %default]",
              dest="markers_all"),
  make_option(c("-n","--duplicatedmarkers"), type="character", default=NULL,
              help="Duplicated markers position [default %default]",
              dest="markers_dup"),
  make_option(c("-s", "--structure"), type="character", default=NULL,
              help="Structure of the input sequence as defined by its compnent sequences [default %default]",
              dest="structure"),
  make_option(c("-q", "--query"), type="character", default=NULL,
              help="Query sequence ID [default %default]",
              dest="queryID")
)

options(error=traceback)
parser <- OptionParser(usage = "%prog -c alignments.coords.txt -s selfmap.coords.txt -o out [options]",option_list=option_list)
params = parse_args(parser)


out_dir = params$out_dir
out_file = params$out_file
# prepare variables
alignments_file=params$coords
alignments_query=params$queryID
alignments_title=paste0(params$queryID," self alignment")

gene_counts_file=params$gene_counts
markers_all_file=params$markers_all
markers_dup_file=params$markers_dup
structure_file=params$structure

category_conversion = c(
"all_genes" = "All genes" ,
"dup_genes" = "Duplicated genes" ,
"all_markers" = "All Markers",
"dup_markers" = "Duplicated Markers"
)

# load data
tb_1 <- read_delim(file = alignments_file,
                         delim = "\t",
                         col_names = c( "tID", "tLen" , "tStart" , "tStop" , "qID" , "qLen" , "qStart" , "qStop" , "identity" , "match" )
                         )

gene_trace <- read_delim(gene_counts_file,
                            delim = "\t",
                            col_names = c( "Chr" , "Start" , "Stop" , "Gene_Id" , "count")
                            )

if (!(empty(gene_trace))) { max_gene_count <- max(gene_trace$count) } else { max_gene_count <- 2 }

all_marker_trace <- read_delim(markers_all_file,
                            delim = "\t",
                            col_names = c( "Chr" , "Start" , "Stop" , "Marker_Id") ,
                            )
dup_markers_trace <- read_delim(markers_dup_file,
                            delim = "\t",
                            col_names = c( "Chr" , "Start" , "Stop" , "Marker_Id") ,
                            )

structure_trace <- read_delim(structure_file,
                            delim = "\t",
                            col_names = c( "Chr" , "Start" , "Stop" , "Seq_Id" , "orientation" ) ,
                            )

coords1 <- tb_1 %>% filter(tID==alignments_query , qID==alignments_query )

refLen <- max(coords1$tLen)
queryLen <- max(coords1$qLen)
self = tibble(tStart = 0 , qStart = 0  , tStop = refLen , qStop = refLen )

# plot
p0 <- ggplot() +
  geom_segment( data = coords1 ,
    aes( x = tStart,
         y = qStart,
         xend = tStop,
         yend = qStop,
         colour = alignments_query ,
        )
  ) +
  geom_segment( data = self ,
    aes( x = tStart,
         y = qStart,
         xend = tStop ,
         yend = qStop ,
         colour = "Self match"
        ),
  ) +
  scale_color_manual(name = "", labels = c("","", ""), breaks = c( alignments_query , "Self match" , "Gene count") , values = c("black", "gray66" , "cornflowerblue") ) +
  scale_x_continuous(limits = c(0, refLen)) +
  scale_y_continuous(limits = c(0, refLen)) +
  labs(title = alignments_title, x = alignments_query, y = alignments_query)

p0 <- p0 + theme(plot.title = element_text(color="black", size=8, face="plain", margin=margin(t=1,r=0,b=1,l=0,"mm")),
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
               plot.margin=unit(c(1,1,1,1),"mm"),
               legend.position="right",
               legend.margin = margin(10,10,10,10) ,
               legend.box.margin = margin(10,10,10,10) ,
               )

p1 <- ggplot(gene_trace) +
  geom_rect(
    aes(
      xmin = Start,
      xmax = Stop,
      ymin = 0 ,
      ymax = count ,
      colour = "Gene count",
      fill = "Gene count",
      group = 1,
      ),
        show.legend = FALSE
       ) +
  scale_color_manual(name = "", breaks = c( alignments_query , "Self match" , "Gene count") , values = c("black", "gray66" , "cornflowerblue") ) +
  scale_fill_manual(name = "", breaks = c( alignments_query , "Self match" , "Gene count") , values = c("black", "gray66" , "cornflowerblue") ) +
  scale_x_continuous(limits = c(0, refLen)) +
  scale_y_continuous(breaks = seq( 0 , max_gene_count, 1)) +
  labs(x = "Sequence position (bp)", y = "Gene counts")

p1 <- p1 +
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
                 plot.margin=unit(c(.7,.7,.7,.7),"mm"),
                 #legend.position="right"
            )



p2 <- ggplot() +
  geom_rect( data=all_marker_trace ,
    aes(
      xmin = Start,
      xmax = Stop,
      ymin = 0 ,
      ymax = 1 ,
      colour = "Marker",
      fill = "Marker",
      group = 1,
      ),
        show.legend = FALSE
       ) +
  scale_color_manual(name = "", breaks = c( "Marker" , "Duplicated marker" ) , values = c("black" ,"red") ) +
  scale_fill_manual(name = "", breaks = c( "Marker" , "Duplicated marker" ) , values = c("black" ,"red") ) +
  scale_x_continuous(limits = c(0, refLen)) +
  scale_y_continuous(breaks = c( 0 , 1)) +
  labs(x = "Sequence position (bp)", y = "Markers")

p2 <- p2 +
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
                 plot.margin=unit(c(.7,.7,.7,.7),"mm"),
                 #legend.position="right"
            )


p4 <- ggplot() +
  geom_rect( data=dup_markers_trace ,
    aes(
      xmin = Start,
      xmax = Stop,
      ymin = 0 ,
      ymax = 1 ,
      colour = "Duplicated marker",
      fill = "Duplicated marker",
      group = 1,
      ),
        show.legend = FALSE
       ) +
  scale_color_manual(name = "", breaks = c( "Marker" , "Duplicated marker" ) , values = c("black" ,"red") ) +
  scale_fill_manual(name = "", breaks = c( "Marker" , "Duplicated marker" ) , values = c("black" ,"red") ) +
  scale_x_continuous(limits = c(0, refLen)) +
  scale_y_continuous(breaks = c( 0 , 1)) +
  labs(x = "Sequence position (bp)", y = "Duplicated<br>Markers")

p4 <- p4 +
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
                 plot.margin=unit(c(.7,.7,.7,.7),"mm"),
                 #legend.position="right"
            )


p3 <- ggplot(structure_trace)

if (!(empty(structure_trace))) {
p3 <- p3 +  geom_rect(
    aes(
      xmin = Start,
      xmax = Stop,
      ymin = 0 ,
      ymax = 1 ,
      #colour = Seq_Id,
      fill = Seq_Id,
      group = 3,
      text = sprintf(
        'Component sequence ID: %s<br>Component region: %s-%s<br>Orientation: %s',
        Seq_Id, Start , Stop, orientation
          )
      ),
    )
} else {
p3 <- p3 +  geom_rect(
    aes(
      xmin = 0,
      xmax = 0,
      ymin = 0 ,
      ymax = 1 ,
      ),
    )
}

p3 <- p3 + labs(colour = "" , fill="") +
  scale_y_continuous(breaks = c( 0 ,1)) +
  scale_x_continuous(limits = c(0, refLen)) +
  labs(x = "Sequence position (bp)", y = "Composition")

p3 <- p3 +
          theme(plot.title = element_text(color="black", size=10, face="plain", margin=margin(t=1,r=0,b=1,l=0,"mm")),
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
                 plot.margin=unit(c(.7,.7,.7,.7),"mm"),
               legend.position="bottom",
               legend.margin = margin(10,10,10,10) ,
               legend.box.margin = margin(10,10,10,10) ,
            )

suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
pdf(NULL)
gt0 <- ggplotGrob(p0)
gt1 <- ggplotGrob(p1)
gt2 <- ggplotGrob(p2)
gt3 <- ggplotGrob(p3)
gt4 <- ggplotGrob(p4)

gt0_widths <- gt0$widths
gt1_widths <- gt1$widths
gt2_widths <- gt2$widths
gt3_widths <- gt3$widths
gt4_widths <- gt4$widths

maxWidth <- unit.pmax(gt0_widths, gt1_widths, gt2_widths, gt3_widths, gt4_widths)

gt0$widths <- maxWidth
gt1$widths <- maxWidth
gt2$widths <- maxWidth
gt3$widths <- maxWidth
gt4$widths <- maxWidth

layout <- rbind(1,2,3,4,5)

g <- arrangeGrob( gt0, gt3, gt1 , gt2 , gt4 , layout_matrix=layout , heights=c(.45, .15, .2 , .1 , .1) )
ggsave(file=paste0(out_dir,"/",out_file, ".pdf"), g , height = 15 , width = 10 )
ggsave(file=paste0(out_dir,"/",out_file, ".png"), g , height = 15 , width = 10 )