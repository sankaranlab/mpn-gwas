library(data.table)
library(tidyverse)
library(BuenColors)
library(plotly)

# Read in full sumstats
sumstats <- "file path to MPN summary statistics"
filtered <- fread(sumstats)

# GFI1B
# Read in LD to sentinel information for lead variant
ld_file <- "file path for LD matrix"
ld_to_lead <- fread(ld_file) %>% 
  mutate(nonsentinel_RSID = ifelse(RSID1=="rs524137",RSID2,RSID1)) %>%
  dplyr::select(nonsentinel_RSID,R2) %>% add_row(., nonsentinel_RSID = "rs524137",R2=1)

chr <- 9
center <- 135879542
width <- 25000
range <- c(center-width,center+width) # the range for rs524137 chr9:135,861,742-135,884,241

lead_variant_pos <- c(135879542)

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

cowplot::ggsave2(p1, file="../output/locus_plots/gfi1b_rs524137_locus_r4.pdf",height = 0.5,width=2)


# CHEK2
# Read in LD to sentinel information for lead variant
ld_file <- "file path for LD matrix"

ld_to_lead <- fread(ld_file) %>% 
  mutate(nonsentinel_RSID = ifelse(RSID1=="rs17879961",RSID2,RSID1)) %>%
  dplyr::select(nonsentinel_RSID,R2) %>% add_row(., nonsentinel_RSID = "rs17879961",R2=1)

chr <- 22
center <- 29121087
width <- 25000
range <- c(center-width,center+width) # the range for rs524137 chr9:135,861,742-135,884,241
lead_variant_pos <- c(29121087)

subset <- sumstats %>% filter(CHR == chr, POS > range[1], POS < range[2])
subset <- merge(subset,ld_to_lead,by.x = "RSID",by.y="nonsentinel_RSID")
subset %>% arrange(pvalue) %>% head()

p1 <- ggplot(subset,aes(x=POS,y=-log10(pvalue),label=MarkerName, label=RSID)) + 
  geom_point(aes(fill=R2),shape=21,size=1)+
  scale_fill_gradientn(colors = jdb_palette("solar_extra")[-c(1,2)],name="R2") +
  guides(fill=guide_colorbar(title.vjust=0.75))+
  labs(x="",y="-log10(p)")+
  geom_vline(xintercept = lead_variant_pos[1], color="red")+
  pretty_plot(fontsize=6) + L_border() +
  scale_y_continuous(expand=c(0,0)) + 
  scale_x_continuous(expand=c(0,0)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position ="none",
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.margin = margin(1, 1, 1, 1),
        legend.direction = "horizontal",
        legend.title = element_text(face="bold")) 

p1

cowplot::ggsave2(p1, file="../output/locus_plots/CHEK2_locus_r4.pdf",height = 1.0,width=2)
