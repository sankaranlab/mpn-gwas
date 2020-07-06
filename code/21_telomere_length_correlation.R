library(tidyverse)
library(data.table)
library(GenomicRanges)
library(BuenColors)
library(qvalue)
library(cowplot)
library(plotly)

sumstats <- "file path to merged MPN and telomere length summary statistics"
merged <- fread(sumstats)

# Subset to locus of interest
finemap <- fread("../data/abf_finemap/MPN_arraycovar_meta_finngen_r4_abf_cojo_PP0.001_annotated.bed")
sentinels <- finemap %>% filter(sentinel=="yes") 
sentinels %>% filter(seqnames=="chr5")
region_range <- 50000
coordinates <- c(1285974-region_range,1285974+region_range)
toplot <- merged %>% filter(CHR==5, POS > coordinates[1], POS < coordinates[2]) %>%
  mutate(sig = ifelse(mpn_pval < 5e-8 & telo_pval < 5e-8, "significant","not significant"),
         z_telo = -beta/se, z_mpn = -Effect/StdErr)
table(toplot$sig)

# Check 2 lead variants
toplot %>% filter(RSID == "rs7705526")
toplot %>% filter(RSID == "rs2853677")

toplot %>% arrange(telo_pval) %>% head()

p2 <- ggplot(toplot,aes(x=z_telo,y=z_mpn)) +
  geom_point(aes(text=RSID,color=sig),size=0.35) +
  geom_smooth(method='lm',color="black",se=FALSE,linetype="dashed",size=0.35) +
  scale_color_manual(values=jdb_palette("flame_light")[c(2,7)]) +
  scale_y_continuous(expand = c(0.01, 0))+
  scale_x_continuous(expand = c(0.01, 0))+
  labs(x="Telomere Length z-score",y = "MPN z-score") + 
  pretty_plot(fontsize = 6) + L_border() +
  theme(plot.title = element_text(hjust = 0.5),legend.position="none") 

cowplot::ggsave2(p2, file="../output/telomere_length/TERT_correlation_r4.pdf",
                 height=3.3,width=3.3,units="cm")
cor.test(toplot$z_telo,toplot$z_mpn,method = "pearson")
