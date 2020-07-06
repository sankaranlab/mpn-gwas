library(data.table)
library(tidyverse)
library(BuenColors)
library(cowplot)

#GEMM" stands for granulocyte, erythrocyte, monocyte, megakaryocyte
cfu_color_maps <- ejc_color_maps[c("Ery","GMP-C","Mono","GMP-A","CMP")]
names(cfu_color_maps) <- c("BFU-E","CFU-G","CFU-M","CFU-GM","CFU-GEMM")

# Read in CFU results
cfu2 <- readxl::read_xlsx("../../data/experimental_data/CFU_assays/Exp2_Day14_primary_CFU-C.xlsx") %>% 
  filter(Sample != "Mock") %>% mutate(Replicate = paste0("exp2_",Replicate))
cfu3 <- readxl::read_xlsx("../../data/experimental_data/CFU_assays/Exp3_CFU-C.xlsx") %>%
  filter(Sample != "Mock") %>% mutate(Replicate = paste0("exp3_",Replicate))

# Combine data from different experiments
cfu <- bind_rows(cfu2,cfu3)

melted <- cfu %>% pivot_longer(-c(Sample,Plating,Replicate),names_to="type",values_to="count")

combined <- melted %>% group_by(Sample,Plating,type) %>% summarise(mean=mean(count),se=sd(count)/sqrt(n()))
combined$Sample <- factor(combined$Sample,levels=c("NT","CDS_g1","CDS_g2"))
combined$type <- factor(combined$type,levels=rev(unique(combined$type)))

errors <- combined %>% group_by(Sample,Plating) %>% mutate(pos = cumsum(mean),  upper = pos + se,  lower = pos - se)
totals <- melted %>% group_by(Sample,Plating,Replicate) %>% summarise(total = sum(count))

p1 <- ggplot(combined,aes(x=Sample,y=mean)) +
  geom_bar(aes(fill = type),color="black", stat="identity") +
  geom_errorbar(data=errors,aes(ymax=upper,  ymin=lower), width=0.15) +
  geom_jitter(data=totals,aes(x=Sample,y=total),color="black",position = position_jitter(width = .2),size=0.8) +
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(1,3,5,7,9)]) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x="",y="No. of colonies") +
  facet_wrap(~Plating) +
  theme(legend.title=element_blank(),axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand=c(0,0))

p1
cowplot::ggsave2(p1, file="../../output/experimental_plots/cfu_assays/gfi1b_cfu_combined.pdf",height=2,width=3)
