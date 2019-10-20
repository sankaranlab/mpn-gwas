library(data.table)
library(tidyverse)
library(BuenColors)
library(cowplot)

# Read in IR results
ir <- readxl::read_xlsx("../../data/experimental_data/CHEK2_data/CHEK2i _HSPC_subpopulations_DO_D4.xlsx")

ir$group <- factor(ir$group,levels=c("DMSO","CHEK2_inh"))
ir$celltype <- factor(ir$celltype,levels=c("CMP","GMP","MEP","HSPC"))

p1 <- ggplot(ir,aes(x=celltype,y=count,fill=group))+ 
  geom_bar(position="dodge", stat="summary",fun.y="mean", color = "black") +
  geom_errorbar(stat = 'summary', fun.data="mean_se",
                width=0.3, position=position_dodge(.9)) +
  geom_jitter(position=position_jitterdodge(dodge.width = 0.9,jitter.width = 0.4),size=0.8) + 
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(1,8)]) +
  labs(x="",y="IR-induced cell death (%)") +
  pretty_plot(fontsize=8) + L_border() +
  scale_y_continuous(expand = c(0, 0.05)) +
  theme(legend.position="none")

cowplot::ggsave2(p1, file="../../output/experimental_plots/chek2/chek2_cycling_IR.pdf",height=1.7,width=2)

# Check for significance
sapply(unique(ir$celltype),function(index){
  print(index)
  
  subset_a <- ir %>% filter(celltype == index,group == "DMSO") %>% .$count
  subset_b <- ir %>% filter(celltype == index,group == "CHEK2_inh") %>% .$count
  
  if (mean(subset_a) == mean(subset_b)){
    return("not significant")
  } else if (length(subset_a) < 2 | length(subset_b) < 2) {
    return("not enough data")
  } else{
    res <- t.test(subset_a,subset_b,alternative=c("two.sided"),paired = TRUE,
                       var.equal = FALSE)
    return(res$p.value)
  } 
})
