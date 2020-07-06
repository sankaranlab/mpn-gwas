library(data.table)
library(tidyverse)
library(BuenColors)
library(ggsignif)

#GEMM" stands for granulocyte, erythrocyte, monocyte, megakaryocyte
cfu_color_maps <- ejc_color_maps[c("Ery","GMP-C","Mono","GMP-A","CMP")]
names(cfu_color_maps) <- c("BFU-E","CFU-G","CFU-M","CFU-GM","CFU-GEMM")

# Read in CFU results
cfu4 <- readxl::read_xlsx("../../data/experimental_data/CFU_assays/Small_enhancer_EXP3_CFU_enh1.xlsx") 

melted <- cfu4 %>% pivot_longer(-c(Sample,Plating,Replicate),names_to="type",values_to="count")

combined <- melted %>% group_by(Sample,Plating,type) %>% summarise(mean=mean(count),se=sd(count)/sqrt(n()))

combined$type <- factor(combined$type,levels=rev(unique(combined$type)))

errors <- combined %>% group_by(Sample,Plating) %>% mutate(pos = cumsum(mean),  upper = pos + se,  lower = pos - se)
totals <- melted %>% group_by(Sample,Plating,Replicate) %>% summarise(total = sum(count))

p2 <- ggplot(combined,aes(x=Sample,y=mean)) +
  geom_bar(aes(fill = type),color="black", stat="identity") +
  geom_errorbar(data=errors,aes(ymax=upper,  ymin=lower), width=0.15) +
  #geom_jitter(data=totals,aes(x=Sample,y=total),color="black",position = position_jitter(width = .2),size=0.8) +
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(1,3,5,7,9)]) +
  pretty_plot(fontsize = 5) + L_border() +
  labs(x="",y="No. of colonies") +
  facet_wrap(~Plating) + 
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(legend.title=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 100),expand=c(0,0)) +
  theme(legend.position = "none")

p2

cowplot::ggsave2(p2, file="../../output/experimental_plots/cfu_assays/gfi1b_cfu_combined_small_enhancer_enh1.pdf",
                 height=3.5,width=3,units="cm")
