library(tidyverse)
library(data.table)
library(cowplot)
library(BuenColors)

gfi1b <- fread("../../data/experimental_data/GFI1B_enhancer_assays/GFI1B_super-enhancer_deletion_genotype.txt")

# Calculate mean and sem of each set of replicates
summaries <- gfi1b %>% group_by(Sample,Outcome) %>% summarise(mean = mean(Percentage),se=sd(Percentage)/sqrt(n()))

errors <- summaries %>% group_by(Sample) %>% mutate(pos = cumsum(mean),  upper = pos + se,  lower = pos - se)
errors$Outcome <- factor(errors$Outcome, levels=c("Uncut","Inversion",'Deletion'))

# merge(gfi1b,errors[,c("Sample","Outcome","pos")],by=c("Sample","Outcome"))
# points <- gfi1b %>% group_by(Sample) %>% mutate(pos = cumsum(Percentage))

p1 <- ggplot(errors,aes(x=Sample,y=mean)) +
  geom_bar(aes(fill = Outcome),color="black", stat="identity",width= 0.7) +
  geom_errorbar(data=errors,aes(ymax=upper,  ymin=lower), width=0.15) +
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(1,3,5,7,9)]) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x="",y="Editing Outcome (%)") +
  theme(legend.title=element_blank(),legend.position="none") +
  scale_y_continuous(expand=c(0,0)) 
p1
cowplot::ggsave2(p1, file="../../output/experimental_plots/gfi1b_enhancer_assays/gfi1b_superenhancer_deletion_genotypes.pdf",
                 height=1.5,width=1.5)
