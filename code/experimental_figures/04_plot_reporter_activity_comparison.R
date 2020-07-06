library(tidyverse)
library(data.table)
library(BuenColors)
library(ggsignif)
library(cowplot)

### reporter assay result comparison

ra_comp <- fread("../../data/experimental_data/reporter_assays/ra_comparison.txt",fill=TRUE)
ra_comp$Sample <- factor(ra_comp$Sample, levels=c("MinP","rs1633768","rs524137"))

ra_comp_p <- ggplot(ra_comp,aes(x=Sample,y=Number,fill=Sample)) +
  geom_bar(position=position_identity(), stat="summary",fun.y="mean", color = "black",width = 0.7) +
  geom_errorbar(stat = 'summary', fun.data="mean_se",
                width=0.15, position=position_dodge(.9),) +
  geom_jitter(size=0.1,color="black",position = position_jitter(width = .1)) +
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(1,8,9)]) +
  labs(x="",y="Relative Luciferase Activity") +
  pretty_plot(fontsize=5) + L_border() + 
  theme(legend.title=element_blank(),legend.key.size = unit(0.5,"line"),
        axis.text.x = element_text(angle = 45,hjust=1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 53),expand=c(0,0)) +
  geom_signif(comparisons = list(c("MinP","rs1633768"),c("MinP","rs524137"),c("rs524137","rs1633768")),
              test = "t.test",y_position = c(10,44,49),textsize=1,size = 0.3,
              annotations = c("ns","**","**"))

ra_comp_p

cowplot::ggsave2(ra_comp_p, file="../../output/experimental_plots/luciferase_assays/reporter_plots/reporter_assay_comparison.pdf",
                 height=2.8,width=1.6,units="cm")
### Double check p values
ra_comp <- as.data.frame(ra_comp)
MinP <- ra_comp %>% filter(Sample == "MinP")
rs1633768 <- ra_comp %>% filter(Sample == "rs1633768")
rs524137 <- ra_comp %>% filter(Sample == "rs524137")


t.test(MinP$Number,rs524137$Number,alternative = "two.sided", var.equal = FALSE)
t.test(rs1633768$Number,rs524137$Number,alternative = "two.sided", var.equal = FALSE)
t.test(rs1633768$Number,MinP$Number,alternative = "two.sided", var.equal = FALSE)
