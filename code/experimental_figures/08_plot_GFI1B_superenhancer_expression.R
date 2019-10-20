library(tidyverse)
library(data.table)
library(cowplot)
library(BuenColors)

gfi1b <- readxl::read_xlsx("../../data/experimental_data/GFI1B_enhancer_assays/gfi1b_SE_D_repeatEXP_gene_expression_day7.xlsx")
gfi1b[,c("group","experiment","replicate")] <- str_split_fixed(gfi1b$Sample, "_",3)

gfi1b <- gfi1b  %>% group_by(group,experiment) %>% summarise(mean=mean(FC),se=sd(FC)/sqrt(n()))
gfi1b$group <- factor(gfi1b$group, levels=c("AAVS","G47","G48"))
gfi1b <- gfi1b %>% filter(group %in% c("AAVS","G48"))

p1 <- ggplot(gfi1b,aes(x=group,y=mean,fill=group))+ 
  geom_bar(position=position_identity(), stat="summary",fun.y="mean", color = "black",width = 0.7) +
  geom_errorbar(stat = 'summary', fun.data="mean_se",
                width=0.3, position=position_dodge(.9)) +
  geom_jitter(color="black",position = position_jitter(width = .2),size=0.8) +
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(1,8)]) +
  labs(x="",y="Relative GFI1B expression") +
  pretty_plot(fontsize=8) + L_border() + 
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.position = "none") 
p1

cowplot::ggsave2(p1,file="../../output/experimental_plots/gfi1b_enhancer_assays/gfi1b_superenhancer_deletion_expression.pdf",
                height=1.5,width=1.5)
