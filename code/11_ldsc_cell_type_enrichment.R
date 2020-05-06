library(tidyverse)
library(data.table)
library(cowplot)
library(BuenColors)
library(qvalue)

# Read in ldsc cts scores
cts <- fread("../data/ldsc/MPN_meta_finngen_r4_ukid_heme_1000GP3_UK10K.cell_type_results.txt")
cts$qvalue <- p.adjust(cts$Coefficient_P_value,method = "fdr")

setorder(cts,Coefficient_P_value)
cts$Name <- factor(cts$Name,levels = cts$Name)

p1 <- ggplot(cts,aes(x=Name,y=-log10(Coefficient_P_value)))+
  geom_bar(width = 1, aes(fill = Name), colour="black",stat = "identity") +
  scale_fill_manual(values = ejc_color_maps) +
  # geom_hline(yintercept = -log10(0.05/nrow(cts)),linetype="dashed") +
  pretty_plot(fontsize=8) + L_border() +
  labs(y="-log10(pvalue)",x="") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))+scale_y_continuous(expand = c(0, 0))
p1

cowplot::ggsave2(p1,file="../output/ldsc/r4_ldsc_cts_heme_atac.pdf",
                width=3,height=2)

