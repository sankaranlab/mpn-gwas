library(data.table)
library(tidyverse)
library(BuenColors)
library(cowplot)

# Read in CFU results
cfu <- readxl::read_xlsx("../../data/experimental_data/CHEK2_data/chek2_inhib_on_CD34+_CFCs.xlsx")

combined <- cfu %>% group_by(Sample,type) %>% summarise(mean=mean(count),se=sd(count)/sqrt(n()))
combined$Sample <- factor(combined$Sample,levels=c("DMSO","CHEK2_inh"))
combined$type <- factor(combined$type,levels=rev(unique(combined$type)))

errors <- combined %>% group_by(Sample) %>% mutate(pos = cumsum(mean),  upper = pos + se,  lower = pos - se)
totals <- cfu %>% group_by(Sample,Experiment) %>% summarise(total = sum(count))

# Total colonies
p1 <- ggplot(combined,aes(x=Sample,y=mean)) +
  geom_bar(aes(fill = type),color="black", stat="identity") +
  geom_errorbar(data=errors,aes(ymax=upper,  ymin=lower), width=0.15) +
  geom_jitter(data=totals,aes(x=Sample,y=total),color="black",position = position_jitter(width = .2),size=0.8) +
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(1,3,5,7,9)]) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x="",y="No. of colonies") +
  theme(legend.position = "none") +
  scale_y_continuous(expand=c(0,0))
p1

# Percents
combined  <- cfu %>% group_by(Sample,Experiment) %>% mutate(percent = count / sum(count)) %>% ungroup() %>%
  group_by(Sample, type) %>%  summarise(mean_percents = mean(percent),se = sd(percent)/n()) %>% ungroup()
combined$Sample <- factor(combined$Sample,levels=c("DMSO","CHEK2_inh"))
combined$type <- factor(combined$type,levels=rev(unique(combined$type)))

errors <- combined %>% group_by(Sample) %>% mutate(pos = cumsum(mean_percents),  upper = pos + se,  lower = pos - se)

p2 <- ggplot(combined,aes(x=Sample,y=mean_percents)) +
  geom_bar(aes(fill = type),color="black", stat="identity") +
  geom_errorbar(data=errors,aes(ymax=upper,  ymin=lower), width=0.15) +
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(1,3,5,7,9)]) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x="",y="Proportion of colonies (%)") +
  geom_text(aes(label = paste0(100*round(mean_percents,2),"%")), 
            position = position_stack(vjust = 0.5), size = 3) +
  theme(legend.position = "none") +
  scale_y_continuous(expand=c(0,0))
p2

cowplot::ggsave2(p1, file="../../output/experimental_plots/chek2/chek2_cfu_counts.pdf",height=2,width=2)
cowplot::ggsave2(p2, file="../../output/experimental_plots/chek2/chek2_cfu_percents.pdf",height=2,width=2)
