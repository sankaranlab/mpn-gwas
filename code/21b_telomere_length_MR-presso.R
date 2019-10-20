library(tidyverse)
library(data.table)
library(BuenColors)
library(cowplot)
library(TwoSampleMR)
library(MRPRESSO)

# Exposure data
telo_sumstats <- as.data.frame(fread("zcat < /Volumes/broad_sankaranlab/telomere_length_GWAS/telomere_length_forMR.txt.gz"))
telo_dat <- format_data(telo_sumstats,type="exposure")

# LD clump
telo_dat <- clump_data(telo_dat,clump_r2=0.001)

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

# Perform MR-presso
# Run MR-PRESSO global method
mr_presso(data = dat,
          BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
          SdOutcome = "se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
          NbDistribution = 1000,  SignifThreshold = 0.05)