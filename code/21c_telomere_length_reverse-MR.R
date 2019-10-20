library(tidyverse)
library(data.table)
library(BuenColors)
library(cowplot)
library(TwoSampleMR)
library(MRPRESSO)

# Exposure data for all SNPs p < 1e-5
mpn_sumstats <- fread("/Volumes/broad_sankaranlab/MPN_GWAS/metaGWAS/Hinds_UKBB_finngen_metal_sumstats_risk_oriented.p1e-5.txt.gz")
exposure_dat <- format_data(mpn_sumstats,type="exposure",
                            snp_col = "RSID",
                            beta_col = "Effect",
                            se_col = "StdErr",
                            effect_allele_col = "risk",
                            other_allele_col = "nonrisk",
                            eaf_col = "RAF",
                            pval_col = "pvalue")

# Outcome data
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/Volumes/broad_sankaranlab/telomere_length_GWAS/ENGAGE_telo_overall_finalrelease.high_sample_size.rsids.txt.gz",
  sep = "\t",
  snp_col = "snp",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effectall",
  other_allele_col = "noneffectall",
  eaf_col = "eaf",
  pval_col = "p_gc",
  samplesize_col = "n"
) %>% mutate(outcome = "telomere length")

# Harmonize data
dat <- harmonise_data(
  exposure_dat = exposure_dat, 
  outcome_dat = outcome_dat
)

# LD clump
dat_clumped <- clump_data(dat,clump_r2=0.001)

# Perform MR
methods_to_use <- c("mr_egger_regression","mr_ivw","mr_weighted_median")
res <- mr(dat_clumped, method_list=methods_to_use)
res

mr_presso(data = dat_clumped,
          BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
          SdOutcome = "se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
          NbDistribution = 1000,  SignifThreshold = 0.05)

# Plots
p1 <- mr_scatter_plot(res, dat)
p1[[1]]

res_single <- mr_singlesnp(dat,single_method="mr_meta_fixed",
                           all_method=methods_to_use)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat,method = TwoSampleMR::mr_ivw)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]
