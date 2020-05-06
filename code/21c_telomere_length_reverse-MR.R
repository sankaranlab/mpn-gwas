library(tidyverse)
library(data.table)
library(BuenColors)
library(cowplot)
library(TwoSampleMR)
library(MRPRESSO)

# Exposure data for all SNPs p < 1e-5
exposure_sumstats <- "file path for exposure sum stats"
mpn_sumstats <- fread(exposure_sumstats)
exposure_dat <- format_data(mpn_sumstats,type="exposure",
                            snp_col = "RSID",
                            beta_col = "Effect",
                            se_col = "StdErr",
                            effect_allele_col = "risk",
                            other_allele_col = "nonrisk",
                            eaf_col = "RAF",
                            pval_col = "pvalue")

# Outcome data
outcome_sumstats <- "file path for outcome sum stats"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outcome_sumstats,
  sep = "\t",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "ea",
  other_allele_col = "oa",
  eaf_col = "eaf",
  pval_col = "pvalue",
  samplesize_col = "n_samples"
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

