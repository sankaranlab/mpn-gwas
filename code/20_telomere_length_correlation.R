library(tidyverse)
library(data.table)
library(GenomicRanges)
library(BuenColors)
library(qvalue)
library(cowplot)

# Plot loci
make_plots <- function(toplot,gene){
  # Plot correlation of pvalues
  p1 <- ggplot(toplot,aes(x=-log10(p_gc),y=-log10(pvalue))) +
    geom_point(aes(color=sig)) +
    scale_color_manual(values=jdb_palette("flame_light")[c(2,7)]) +
    scale_y_continuous(expand = c(0.01, 0))+
    scale_x_continuous(expand = c(0.01, 0))+
    labs(x="Telomere Length -log10(p)",y = "MPN -log10(p)") + 
    pretty_plot() + L_border() +
    theme(plot.title = element_text(hjust = 0.5),legend.position="none") 
  p2 <- ggplot(toplot,aes(x=z_telo,y=z_mpn)) +
    geom_point(aes(color=sig)) +
    geom_smooth(method='lm',color="black",se=FALSE,linetype="dashed") +
    scale_color_manual(values=jdb_palette("flame_light")[c(2,7)]) +
    scale_y_continuous(expand = c(0.01, 0))+
    scale_x_continuous(expand = c(0.01, 0))+
    labs(x="Telomere Length z-score",y = "MPN z-score") + 
    pretty_plot() + L_border() +
    theme(plot.title = element_text(hjust = 0.5),legend.position="none") 
  
  plot_grid(p1, p2, labels = c('a', 'b'), ncol = 2)
 
}

# TERT 250kb region
toplot <- fread("/Volumes/broad_sankaranlab/MPN_GWAS/metaGWAS/TERT_locus.txt") %>%
  mutate(sig = ifelse(pvalue & p_gc < 5e-8, "significant","not significant"),
         z_telo = -beta/se, z_mpn = -Effect/StdErr)

# Check 2 lead variants
toplot %>% filter(RSID == "rs7705526")
toplot %>% filter(RSID == "rs2853677")

final_plot <- make_plots(toplot,gene="TERT")
cowplot::ggsave2(final_plot, file="../output/telomere_length/TERT_correlation.pdf",
                height=2.5,width=5)
cor.test(toplot$z_telo,toplot$z_mpn,method = "pearson")

gene = "TERT"
p2 <- ggplot(toplot,aes(x=z_telo,y=z_mpn)) +
  geom_point(aes(color=sig)) +
  geom_smooth(method='lm',color="black",se=FALSE,linetype="dashed") +
  scale_color_manual(values=jdb_palette("flame_light")[c(2,7)]) +
  scale_y_continuous(expand = c(0.01, 0))+
  scale_x_continuous(expand = c(0.01, 0))+
  labs(x="Telomere Length z-score",y = "MPN z-score") + 
  pretty_plot(fontsize=8) + L_border() +
  theme(plot.title = element_text(hjust = 0.5),legend.position="none") 
cowplot::ggsave2(p2, file="../output/telomere_length/TERT_zscores.pdf",
                height=2,width=2)

