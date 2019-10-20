library(tidyverse)
library(data.table)
library(BuenColors)
library(cowplot)
library(TwoSampleMR)
library(MendelianRandomization)

# Exposure data (pval < 1e-5)
telo_sumstats <- as.data.frame(fread("zcat < /Volumes/broad_sankaranlab/telomere_length_GWAS/telomere_length_forMR.txt.gz"))
telo_dat <- format_data(telo_sumstats,type="exposure")

# Outcome data
outcome_dat <- read_outcome_data(
  snps = telo_dat$SNP,
  filename = "/Volumes/broad_sankaranlab/MPN_GWAS/metaGWAS/Hinds_UKBB_finngen_metal_sumstats_risk_oriented.txt.gz",
  sep = "\t",
  snp_col = "RSID",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "risk",
  other_allele_col = "nonrisk",
  eaf_col = "RAF",
  pval_col = "pvalue"
) %>% mutate(outcome = "MPN",samplesize.outcome=758103)

# Harmonize data
dat <- harmonise_data(
  exposure_dat = telo_dat, 
  outcome_dat = outcome_dat
)

# LD clump
dat_clumped <- clump_data(dat,clump_r2=0.001)

# Perform MR
mr_method_list()
methods_to_use <- c("mr_egger_regression","mr_ivw","mr_weighted_median")
res <- mr(dat_clumped, method_list=methods_to_use)
res

mr_presso(data = dat_clumped,
          BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
          SdOutcome = "se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
          NbDistribution = 1000,  SignifThreshold = 0.05)

# Plots
p1 <- mr_scatter_plot(res, dat_clumped)
p1 <- p1[[1]] + pretty_plot(fontsize=7) + L_border() + theme(legend.position ="top") +
  geom_hline(yintercept=0,linetype="dashed")
p1

res_single <- mr_singlesnp(dat_clumped,single_method="mr_meta_fixed",
                           all_method=methods_to_use)
p2 <- mr_forest_plot(res_single)
p2 <- p2[[1]]+ pretty_plot(fontsize=7) + L_border() + theme(legend.position = "none")

res_loo <- mr_leaveoneout(dat_clumped,method = TwoSampleMR::mr_ivw)
p3 <- mr_leaveoneout_plot(res_loo)
p3 <- p3[[1]]+ pretty_plot(fontsize=7) + L_border() + theme(legend.position = "none")

# Save plots 
if (TRUE){
  cowplot::ggsave2(p1, file="../output/telomere_length/mendelian_randomization/MR_scatterplot.pdf", 
                  width=2.5, height=2.5)
  cowplot::ggsave2(p2, file="../output/telomere_length/mendelian_randomization/MR_singleSNP.pdf", 
                  width=2.5, height=3)
  cowplot::ggsave2(p3, file="../output/telomere_length/mendelian_randomization/MR_leave_one_out.pdf", 
                  width=2.5, height=3)
}