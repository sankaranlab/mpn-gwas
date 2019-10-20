library(tidyverse)
library(data.table)
library(qvalue)
library(BuenColors)
library(cowplot)

# Read 500K UKBB summary statistics
gencors_MPN <- fread("../data/ldsc/MPN_meta_finngen_ukid_vs_bloodtraits_ukbb_500k_constrained.gc.summary.txt") %>%
  dplyr::select(p2,rg,se,z,p)


gencors_MPN <- gencors_MPN[complete.cases(gencors_MPN),]
gencors_MPN$qvalue <- p.adjust(gencors_MPN$p,method="fdr")
gencors_MPN$significance <- ifelse(gencors_MPN$qvalue < 0.05,"FDR","no")
gencors_MPN[gencors_MPN$qvalue>=0.05,]$significance <- ifelse(gencors_MPN[gencors_MPN$qvalue>=0.05,]$p < 0.05,
                                                        "nominal","no")

pval_colormap <- jdb_palette("Zissou")[c(1,3,5)]; 
names(pval_colormap) <- c("no","nominal","FDR")

gencors_MPN$p2 <- gsub("_"," ",gencors_MPN$p2) %>% gsub("\\..*"," ",.)
gencors_MPN$significance <- factor(gencors_MPN$significance,levels=c("no","nominal","FDR"))

mpn_gc <- ggplot(gencors_MPN,aes(x=p2,color=significance)) +
  geom_pointrange(aes(y = rg, ymin = rg-se, ymax = rg+se),size=0.25) +
  coord_flip() +
  scale_color_manual(values=pval_colormap) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x="",y="MPN genetic correlation") +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position="bottom")
mpn_gc

cowplot::ggsave2(mpn_gc, file = "../output/ldsc/genetic_correlations_mpn_constrained.pdf", width = 2.5, height =2.75)

# gencors with other blood cancers ----------------------------------------

gencors_MPN <- fread("../data/ldsc/MPN_meta_finngen_ukid_vs_otherphenos_SAIGE.gc.summary.txt") 
gencors_MPN <- gencors_MPN[complete.cases(gencors_MPN)]
gencors_MPN$qvalue <- qvalue(gencors_MPN$p,lambda=0)$qvalue
gencors_MPN$significance <- ifelse(gencors_MPN$qvalue < 0.05,"yes","no")
gencors_MPN[gencors_MPN$qvalue>=0.05,]$significance <- ifelse(gencors_MPN[gencors_MPN$qvalue>=0.05,]$p < 0.05,
                                                              "nominal","no")
gencors_MPN$p2 <- gsub(".sumstats.gz","",gencors_MPN$p2) %>% gsub("_"," ",.)

pval_colormap <- jdb_palette("Zissou")[c(1,3,5)]; 
names(pval_colormap) <- c("no","nominal","yes")

mpn_gc <- ggplot(gencors_MPN,aes(x=p2,color=significance)) +
  geom_pointrange(aes(y = rg, ymin = rg-se, ymax = rg+se)) +
  coord_flip() +
  scale_color_manual(values=pval_colormap) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x="",y="MPN genetic correlation") +
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position="none")
cowplot::ggsave2(mpn_gc, file = "../output/ldsc/genetic_correlations_mpn_otherphenotypes.pdf", width = 2.5, height = 2)

