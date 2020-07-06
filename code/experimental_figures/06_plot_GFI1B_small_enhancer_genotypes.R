library(tidyverse)
library(data.table)
library(BuenColors)
library(ggsignif)
library(cowplot)

### GFI1B Small Enh Deletion EXP Genotype


geno <- fread("SankaranLab/mpn-GWAS/data/experimental_data/GFI1B_small_enh_deletion/genotype_summary.txt")

# Calculate mean and sem of each set of replicates
summaries <- geno  %>% group_by(Sample,Outcome) %>% summarise(mean = mean(Percentage),se=sd(Percentage)/sqrt(n()))

errors <- summaries %>% group_by(Sample) %>% mutate(pos = cumsum(mean),  upper = pos + se,  lower = pos - se)
errors$Outcome <- factor(errors$Outcome, levels=c("Uncut","Edited"))


p1 <- ggplot(errors,aes(x=Sample,y=mean)) +
  geom_bar(aes(fill = Outcome),color="black", stat="identity",width= 0.7) +
  geom_errorbar(data=errors,aes(ymax=upper,  ymin=lower), width=0.15) +
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(1,3,7)]) +
  pretty_plot(fontsize = 6) + L_border() +
  labs(x="",y="Editing Outcome (%)") +
  theme(legend.title=element_blank(),legend.position="right",legend.key.size = unit(1.0,"line"),
        axis.ticks = element_line(colour = "black", size = 1), axis.text.x = element_text(angle = 45,hjust=1)) +
  scale_y_continuous(expand=c(0,0)) 

p1

cowplot::ggsave2(p1, file="SankaranLab/mpn-GWAS/output/experimental_plots/gfi1b_small_enhancer_deletion/genotype.pdf",
                 height=1.5,width=1.5)
