library(tidyverse)
library(data.table)
library(BuenColors)
library(ggsignif)
library(cowplot)

### reporter assay

exp4 <- fread("../..//data/experimental_data/reporter_assays/EXP4_rs524137.txt",fill=TRUE)

exp4_p <- ggplot(exp4,aes(x=Sample,y=Mean,fill=Sample)) +
  geom_bar(position=position_identity(), stat="summary",fun.y="mean", color = "black",width = 0.7) +
  geom_errorbar(stat = 'summary', fun.data="mean_se",
                width=0.15, position=position_dodge(.9),) +
  geom_jitter(size=0.3,color="black",position = position_jitter(width = .1)) +
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(1,8,9)]) +
  labs(x="",y="GFP MFI") +
  pretty_plot(fontsize=5) + L_border() + 
  theme(axis.ticks = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 45,hjust=1),
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 120000),expand=c(0,0)) +
  # ggtitle("rs524137") +
  # geom_signif(comparisons = list(c("Non-risk","Risk")),
  #             test = "t.test",y_position = 105000,textsize = 2,size = 0.3,
  #             annotations ="***") +
  theme(legend.position="none")

cowplot::ggsave2(exp4_p, file="../../output/experimental_plots/luciferase_assays/reporter_plots/EXP4_rs524137.pdf",
                 height=2.6,width=1.6,units="cm")

### Double check p values
exp4 <- as.data.frame(exp4)
r <- exp4[exp4$Sample=="Risk",]
nr <- exp4[exp4$Sample=="Non-risk",]
t.test(r$Mean,nr$Mean,alternative = "two.sided", var.equal = FALSE)
