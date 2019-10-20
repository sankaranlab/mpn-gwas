library(tidyverse)
library(data.table)
library(cowplot)
library(BuenColors)

# Read in reporter results
gfi1b <- fread("../../data/experimental_data/reporter_assays/GFI1B_reporter_results.txt")
gfi1b$type <- factor(gfi1b$type, levels=c("MinP","nonrisk","risk"))

plot_and_save <- function(reporters,name){
  # Plot
  p1 <- ggplot(reporters,aes(x=type,y=activity,fill=type))+ 
    geom_bar(position=position_identity(), stat="summary",fun.y="mean", color = "black", width = 0.7) +
    geom_errorbar(stat = 'summary', fun.data="mean_se",
                  width=0.3, position=position_dodge(.9)) +
    geom_jitter(color="black",position = position_jitter(width = .2),size=0.8) +
    scale_fill_manual(values = jdb_palette("brewer_spectra")[c(6,1,8)]) +
    labs(x="",y="Relative luciferase activity") +
    pretty_plot(fontsize=8) + L_border() +
    scale_y_continuous(expand = c(0, 0))+
    theme(legend.position = "none") 
  
  cowplot::ggsave2(p1,file=paste0("../../output/experimental_plots/luciferase_assays/reporter_plots/reporter_barplot_",name,".pdf"),
                  height=1.5,width=1.5)
  return(p1)
}

plot_and_save(gfi1b,"GFI1B")
