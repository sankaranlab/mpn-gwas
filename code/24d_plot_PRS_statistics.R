library(data.table)
library(tidyverse)
library(BuenColors)
library(ggrastr)

# Read in PRS table
scalar <- "1e-5"
all_covariates <- fread(paste0("filepath")) %>% 
  mutate(pheno = ifelse(pheno == 0,"Control","Case"))

# Plot PRS distribution
d <- ggplot(all_covariates, aes(PRS)) +
  geom_density() +
  pretty_plot(fontsize = 7)+ L_border() +
  labs(x="Polygenic score",y="Density") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(xintercept = median(all_covariates$PRS), linetype = 2)

cowplot::ggsave2(d, file=paste0("../output/PRS/r4_PRS_",scalar,"_distribution.pdf"),width=1.5,height=1.5)

# Plot cases vs. controls
all_covariates %>% group_by(pheno) %>% summarise(median(percentile_PRS))

c <- ggplot(all_covariates, aes(y=percentile_PRS,x=pheno,color=pheno,fill=pheno)) +
  geom_boxplot(varwidth = F,alpha = 0.6,notch = TRUE) +
  scale_color_manual(values=jdb_palette("brewer_spectra")[c(7,1)]) +
  pretty_plot(fontsize = 7)+
  L_border() +
  labs(x="",y="Polygenic score percentile") +
  theme(legend.position="none")
c

cowplot::ggsave2(c, file=paste0("../output/PRS/r4_PRS_",scalar,"_case_vs_control.pdf"),width=1.5,height=1.5)

  