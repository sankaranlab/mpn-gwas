library(tidyverse)
library(data.table)
library(BuenColors)
library(ggsignif)
library(cowplot)

### GFI1B Small Enh Deletion EXP mRNA

mRNA <- fread("../../data/experimental_data/GFI1B_small_enh_deletion/mRNA_expression_ACTIN_combined.txt",fill=TRUE)

mRNA$Type <- factor(mRNA$Type, levels=c("bulk","sorted"))
mRNA <- mRNA %>% filter(Type %in% c("bulk","sorted"))

mRNA$Sample <- factor(mRNA$Sample,levels=c("AAVS1_bulk","ENH_bulk","AAVS1_sorted","ENH_sorted"))

p1 <- ggplot(mRNA,aes(x=Sample,y=Mean,fill=Type)) +
  geom_bar(position=position_identity(), stat="summary", color = "black",width = 0.7) +
  geom_errorbar(stat = 'summary', fun.data="mean_se",
                width=0.3, position=position_dodge(.9)) +
  geom_jitter(size=0.3,color="black",position = position_jitter(width = .1)) +
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(1,7)]) +
  labs(x="",y="Relative"~italic("GFI1B")~"expression") +
  pretty_plot(fontsize=6) + L_border() + 
  theme(legend.title=element_blank(),legend.position="right",legend.key.size = unit(1,"line"),
        axis.ticks = element_line(colour = "black", size = 1), axis.text.x = element_text(angle = 45,hjust=1)) +
  scale_y_continuous(limits = c(0, 2.1),expand=c(0,0)) +
  geom_signif(comparisons = list( c("AAVS1_bulk","ENH_bulk"),
                                  c("AAVS1_sorted","ENH_sorted")),
              test = "t.test",y_position = c(1.95,1.95),textsize = 2,size = 0.3,
              annotations =NULL) +
  scale_x_discrete(labels = c('AAVS1','ENH','AAVS1','ENH')) +
  theme(legend.position = "none")

p1
cowplot::ggsave2(p1, file="../../output/experimental_plots/gfi1b_small_enhancer_deletion/mRNA.pdf",
                 height=1.5,width=2)

### Double check p values
mRNA <- as.data.frame(mRNA)
Ab <- mRNA[mRNA$Sample=="AAVS1_bulk",]
Eb <- mRNA[mRNA$Sample=="ENH_bulk",]
t.test(Ab$Mean,Eb$Mean,alternative = "two.sided", var.equal = FALSE)
As <- mRNA[mRNA$Sample=="AAVS1_sorted",]
Es <- mRNA[mRNA$Sample=="ENH_sorted",]
t.test(As$Mean,Es$Mean,alternative = "two.sided", var.equal = FALSE)

