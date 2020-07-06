library(tidyverse)
library(data.table)
library(BuenColors)
library(ggsignif)
library(cowplot)

### GFI1B Small Enh Deletion EXP Cell Growth FC
raw_counts <- fread("../../data/experimental_data/GFI1B_small_enh_deletion/growth_number.txt")
  
  
FC <- fread("../../data/experimental_data/GFI1B_small_enh_deletion/growth_fc.txt",fill=TRUE) %>%
  mutate(log2FC = log2(FC), Sample = ifelse(Sample == "Total","All cells",Sample))

FC$Sample <- factor(FC$Sample, levels=c("All cells","CD34","LT-HSC"))
FC <- as.data.frame(FC)

FCp<- ggplot(FC,aes(x=Sample,y=log2FC,fill=Sample)) +
  geom_bar(position=position_identity(), stat="summary",fun.y="mean", color = "black",width = 0.7) +
  geom_errorbar(stat = 'summary', fun.data="mean_se",
                width=0.3, position=position_dodge(.9)) +
  geom_jitter(size=0.2,color="black",position = position_jitter(width = .1)) +
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(1,8,9)]) +
  labs(x="",y="Log2 Fold Change Relative to AAVS1") +
  pretty_plot(fontsize=6) + L_border() + 
  theme(legend.title=element_blank(),legend.position=NULL,legend.key.size = unit(0.8,"line"),
        axis.ticks = element_line(colour = "black", size = 1), 
        axis.text.x = element_text(angle = 45,hjust=1),
        plot.title = element_text(hjust = 0.5,size=4.5,face="bold")) +
  scale_y_continuous(limits = c(-0.1, 3),expand=c(0,0)) +
  # ggtitle("log2(FC) in Cell Numbers") +
  theme(legend.position = "none") +
  geom_signif(comparisons = list(c("All cells","LT-HSC"),c("All cells","CD34"),c("LT-HSC","CD34")),
              test = "t.test",y_position = c(2.8,1.3,2.5),textsize = 1.5,size = 0.3,
              annotations = NULL,test.args=list(alternative = "two.sided", var.equal = FALSE, paired=TRUE)) +
  theme(legend.position = "none")

FCp


cowplot::ggsave2(FCp, file="../../output/experimental_plots/gfi1b_small_enhancer_deletion/growth_fold_change_log2.pdf",
                 height=1.5,width=1.0)

### Double check p values
fc <- as.data.frame(FC)
total <- fc[fc$Sample=="Total",]
hsc <- fc[fc$Sample=="LT-HSC",]
cd34 <- fc[fc$Sample=="CD34",]

t.test(total$log2FC,hsc$log2FC,alternative = "two.sided", var.equal = FALSE,paired = TRUE)
t.test(total$log2FC,cd34$log2FC,alternative = "two.sided", var.equal = FALSE,paired = TRUE)
t.test(cd34$log2FC,hsc$log2FC,alternative = "two.sided", var.equal = FALSE,paired = TRUE)


t.test(total$FC,hsc$FC,alternative = "two.sided", var.equal = FALSE,paired = TRUE)
t.test(total$FC,cd34$FC,alternative = "two.sided", var.equal = FALSE,paired = TRUE)
t.test(cd34$FC,hsc$FC,alternative = "two.sided", var.equal = FALSE,paired = TRUE)
