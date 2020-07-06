library(tidyverse)
library(data.table)
library(cowplot)
library(BuenColors)

# Read in reporter results
chek2 <- readxl::read_xlsx("../../data/experimental_data/CHEK2_data/shCHK2_expansion_raw_data.xlsx") 
chek2$type <- factor(chek2$type, levels=c("shCtr","shCHK2"))

sapply(unique(chek2$day),function(day_index){
  print(day_index)
  
  # Group by experiment
  if (FALSE){
    subset_a <- chek2 %>% filter(day == day_index,type == "shCtr") %>% group_by(experiment) %>% summarise(n_control = mean(count)) 
    subset_b <- chek2 %>% filter(day == day_index,type == "shCHK2") %>% group_by(experiment) %>% summarise(n_shchek2 = mean(count)) 
    
    merged <- merge(subset_a,subset_b,by="experiment")
    merged <- merged[complete.cases(merged),]
    subset_a <- merged$n_control
    subset_b <- merged$n_shchek2
    
    if (nrow(merged) < 1){ return("not enough data")}
  }
  
  subset_a <- chek2 %>% filter(day == day_index,type == "shCtr") %>% .$count
  subset_b <- chek2 %>% filter(day == day_index,type == "shCHK2") %>% .$count
  
  if (mean(subset_a) == mean(subset_b)){
    return("not significant")
  } else if (length(subset_a) < 2 | length(subset_b) < 2) {
    return("not enough data")
  } else{
    res <- wilcox.test(subset_a,subset_b,alternative=c("two.sided"),paired = FALSE,
                       var.equal = FALSE)
    return(res$p.value)
  } 
})

means <- chek2 %>% group_by(day,type) %>% summarise(mean = mean(count),sem = sd(count)/sqrt(n()))

# Plot
require(scales)
p1 <- ggplot(means,aes(x=day,y=mean,color=type))+ 
  geom_point(size=0.8) + geom_line()+
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),
                 position=position_dodge(.9)) +
  scale_color_manual(values = jdb_palette("brewer_spectra")[c(1,8)]) +
  scale_y_log10(limits = c(0.1,100000),breaks=c(0.1,1,10,100,100,1000,10000,100000),labels=comma)+
  labs(x="Days",y="Expansion") +
  pretty_plot(fontsize=5) + L_border() +
  theme(legend.title = element_blank(),legend.position = "none")

cowplot::ggsave2(p1,file="../../output/experimental_plots/chek2/chek2_expansion.pdf",
                height=3.5,width=3,units="cm")
