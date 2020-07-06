library(tidyverse)
library(data.table)
library(BuenColors)
library(ggsignif)
library(cowplot)

### GFI1B Small Enh Deletion EXP Cell Growth FC

Number <- fread("../data/experimental_data/GFI1B_small_enh_deletion/growth_number.txt",fill=TRUE)

Number$Type <- factor(Number$Type, levels=c("Total","CD34","LT-HSC"))
Number <- Number %>% filter(Type %in% c("Total","CD34","LT-HSC"))

Number$Sample <- factor(Number$Sample,levels=c("AAVS1_total","ENH_total","AAVS1_CD34","ENH_CD34",
                                             "AAVS1_LT-HSC","ENH_LT-HSC"))

p1 <- ggplot(Number,aes(x=Sample,y=Number,fill=Type)) +
  geom_bar(position=position_identity(), stat="summary",fun.y="mean", color = "black",width = 0.7) +
  geom_errorbar(stat = 'summary', fun.data="mean_se",
                width=0.3, position=position_dodge(.9)) +
  geom_jitter(size=0.3,color="black",position = position_jitter(width = .1)) +
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(1,5,7)]) +
  labs(x="",y="Cell Numbers") +
  pretty_plot(fontsize=6) + L_border() + 
  theme(legend.title=element_blank(),legend.position="right",legend.key.size = unit(1,"line"),
        axis.ticks = element_line(colour = "black", size = 1), axis.text.x = element_text(angle = 45,hjust=1)) +
  scale_y_continuous(limits = c(0, 4000000),expand=c(0,0)) +
  geom_signif(comparisons = list( c("AAVS1_total","ENH_total"),
                                  c("AAVS1_CD34","ENH_CD34"),
                                  c("AAVS1_LT-HSC","ENH_LT-HSC")),
              test = "t.test",y_position = c(3800000,3500000,500000),textsize = 2,size = 0.3,
              annotations =NULL,test.args=list(alternative = "two.sided", var.equal = FALSE, paired=FALSE)) +
  scale_x_discrete(labels = c('AAVS1','ENH','AAVS1','ENH','AAVS1','ENH')) +
  theme(legend.position = "none") 

p1


cowplot::ggsave2(p1, file="../output/experimental_plots/gfi1b_small_enhancer_deletion/growth_number.pdf",
                 height=1.5,width=1.0)
### Double check p values
number <- as.data.frame(Number)
total <- number[number$Type=="Total",]
hsc <- number[number$Type=="LT-HSC",]
cd34 <- number[number$Type=="CD34",]
t.test(hsc$Number[1:3],hsc$Number[4:9],alternative = "two.sided", var.equal = FALSE)
t.test(cd34$Number[1:3],cd34$Number[4:9],alternative = "two.sided", var.equal = FALSE)
t.test(total$Number[1:3],total$Number[4:9],alternative = "two.sided", var.equal = FALSE)