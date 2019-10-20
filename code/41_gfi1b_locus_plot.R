library(data.table)
library(tidyverse)
library(BuenColors)
library(plotly)

# Read in full sumstats
sumstats <- fread("zcat < /Volumes/broad_sankaranlab/MPN_GWAS/metaGWAS/Hinds_UKBB_finngen_metal_sumstats_risk_oriented.txt.gz")

# Read in LD to sentinel information for lead variant
ld_to_lead <- fread("../data/gfi1b_locus/ukb_region.rs621940.marked.ld") %>% 
  mutate(nonsentinel_RSID = ifelse(RSID1=="rs621940",RSID2,RSID1)) %>%
  dplyr::select(nonsentinel_RSID,R2) %>% add_row(., nonsentinel_RSID = "rs621940",R2=1)

# GFI1B
chr <- 9
center <- 135870130
width <- 50000
range <- c(135848632,135889482) # the range for Figure 3E
# range <- c(135849949,135884754) # the range for Figure S11
lead_variant_pos <- c(135870130)

subset <- sumstats %>% filter(CHR == chr, POS > range[1], POS < range[2])
subset <- merge(subset,ld_to_lead,by.x = "RSID",by.y="nonsentinel_RSID")
subset %>% arrange(pvalue) %>% head()

p1 <- ggplot(subset,aes(x=POS,y=-log10(pvalue),label=MarkerName, label=RSID)) + 
  geom_point(aes(fill=R2),shape=21,size=1)+
  scale_fill_gradientn(colors = jdb_palette("solar_extra")[-c(1,2)],name="R2") +
  guides(fill=guide_colorbar(title.vjust=0.75))+
  labs(x="",y="-log10(p)")+
  geom_hline(yintercept = -log10(5e-8),linetype="dashed")+
  geom_vline(xintercept = lead_variant_pos[1], color="red")+
  pretty_plot(fontsize=6) + L_border() +
  scale_y_continuous(expand=c(0,0)) + 
  scale_x_continuous(expand=c(0,0)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(.05, .9),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.margin = margin(1, 1, 1, 1),
        legend.direction = "horizontal",
        legend.title = element_text(face="bold")) 
  

p1
# ggplotly(p1)

cowplot::ggsave2(p1, file="../output/locus_plots/gfi1b_locus.pdf",height = 0.5,width=2)

