
library(tidyverse)
library(data.table)
library(BuenColors)
library(ggsignif)
library(cowplot)

### GFI1B Small Enh Deletion EXP FACS LT-HSC Number

sorted <- fread("../../data/experimental_data/GFI1B_small_enh_deletion/final_sorted_total_cell_number_combined.txt",fill=TRUE)

p1 <- ggplot(sorted,aes(x=Sample,y=Number,fill=Sample)) +
  geom_bar(position=position_identity(), stat="summary", color = "black",width = 0.7) +
  geom_errorbar(stat = 'summary', fun.data="mean_se",
                width=0.3, position=position_dodge(.9)) +
  geom_jitter(size=0.3,color="black",position = position_jitter(width = .1)) +
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(1,8,9)]) +
  labs(x="",y="No. of Phenotypic LT-HSCs") +
  pretty_plot(fontsize=6) + L_border() + 
  theme(legend.title=element_blank(),legend.position="right",legend.key.size = unit(0.8,"line"),
        axis.ticks = element_line(colour = "black", size = 1), axis.text.x = element_text(angle = 45,hjust=1)) +
  scale_y_continuous(limits = c(0, 90000),expand=c(0,0)) +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45,hjust=1)) +
  #scale_y_continuous(limits = c(0, 2.3),expand=c(0,0)) +
  geom_signif(comparisons = list(c("AAVS1","ENH")),
              test = "t.test",y_position = 81000,textsize =3,size = 0.3,
              annotations = "*") + theme(legend.position = "none")

p1

cowplot::ggsave2(p1, file="../../output/experimental_plots/gfi1b_small_enhancer_deletion/LT-HSC_perc.pdf",
                 height=1.5,width=1.0)

##### Double check p values
sorted <-  as.data.frame(sorted)
A1 <- sorted[sorted$Sample=="AAVS1",]
ENH <- sorted[sorted$Sample=="ENH",]
t.test(A1$Number,ENH$Number,alternative = "two.sided", paired = FALSE, var.equal = FALSE)



