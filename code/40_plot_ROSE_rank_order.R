library(ggplot2)
library(scales)
library(tidyverse)
library(cowplot)
library(data.table)
library(BuenColors)
library(ggrepel)

# Read in MACS peaks file to browse histone peaks (not needed for figure generation)
peaks <- fread("../data/histone_modifications/HSPC_H3K27Ac_peaks.narrowPeak") 
colnames(peaks) <- c("chr","start","end","id","int_score",".","FC","log10p","log10q","summit_position_to_start")

# Read in ROSE
ROSE <- readxl::read_xlsx("../data/histone_modifications/rose_o/HSPC_H3K27Ac_peaks_AllEnhancers_table_gfi1b.xlsx")  %>%
  dplyr::rename(signal = "HSPC_H3K27Ac.no_duplicates.bam")

super.threshold <- 9730.1403 # extracted from ROSE output txt file

callout <- "GFI1B (241)"
ROSE$callout[241] <- callout

p1 <- ggplot(ROSE, aes(x=enhancerRank, y=signal)) + 
  geom_hline(yintercept = super.threshold,colour="grey", linetype = "longdash") +
  geom_point(data=ROSE,color = "#5E4FA2",size=0.6)  + # find pretty color: jdb_palette("brewer_spectra")[c(1,8)]
  geom_point(data=ROSE[which(ROSE$callout == callout),],aes(x=enhancerRank, y=signal),color='#DC494C',size=1) +
  labs(x = "Enhancer Ranking") +labs(y = "H3K27ac Signal") + 
  pretty_plot(fontsize = 8) + L_border() +
  geom_text_repel(data=ROSE,label=ROSE$callout,color = '#DC494C',hjust=1.5,vjust=0.5,aes(fontface=3))+
  theme(aspect.ratio=1) +
  scale_x_reverse() + scale_y_continuous(labels = comma) 
#p1

cowplot::ggsave2(p1,file="../output/rose/rose_rank_order.pdf",height=3,width=3)

