library(tidyverse)
library(data.table)
library(BuenColors)
library(cowplot)
library(TwoSampleMR)
library(MendelianRandomization)
library(MRPRESSO)

# Exposure data (pval < 1e-5)
sumstats <- "file path to telomere length summary statistics"
telo_sumstats <- as.data.frame(fread(sumstats)) %>%
  mutate(Phenotype = "telomere_length")

telo_dat <- format_data(telo_sumstats,type="exposure",
                        snp_col = "rsid",
                        beta_col = "beta", se_col = "se",
                        eaf_col = "eaf",
                        effect_allele_col = "ea",
                        other_allele_col = "oa",
                        pval_col = "pvalue",
                        samplesize_col = "n_samples") 

# Outcome data
outcome_sumstats <- "file path to outcome summary statistics"
outcome_dat <- read_outcome_data(
  snps = telo_dat$SNP,
  filename = outcome_sumstats,
  sep = "\t",
  snp_col = "RSID",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "risk",
  other_allele_col = "nonrisk",
  eaf_col = "RAF",
  pval_col = "pvalue"
) %>% mutate(outcome = "MPN",samplesize.outcome=838503)

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
mr_scatter_plot <- function(mr_results, dat)
{
  # dat <- subset(dat, paste(id.outcome, id.exposure) %in% paste(mr_results$id.outcome, mr_results$id.exposure))
  requireNamespace("ggplot2", quietly=TRUE)
  requireNamespace("plyr", quietly=TRUE)
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(nrow(d) < 2 | sum(d$mr_keep) == 0)
    {
      return(blank_plot("Insufficient number of SNPs"))
    }
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
    mrres$a <- 0
    if("MR Egger" %in% mrres$method)
    {
      temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }
    
    if("MR Egger (bootstrap)" %in% mrres$method)
    {
      temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
    }
    
    ggplot2::ggplot(data=d, ggplot2::aes(x=beta.exposure, y=beta.outcome)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0,size=0.35) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0,size=0.35) +
      ggplot2::geom_point(ggplot2::aes(text=paste("SNP:", SNP)),size=0.35) +
      ggplot2::geom_abline(data=mrres, ggplot2::aes(intercept=a, slope=b, colour=method),size=0.35, show.legend=TRUE) +
      ggplot2::scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
      ggplot2::labs(colour="MR Test", x=paste("SNP effect on", d$exposure[1]), y=paste("SNP effect on", d$outcome[1])) +
      ggplot2::theme(legend.position="top", legend.direction="vertical") +
      ggplot2::guides(colour=ggplot2::guide_legend(ncol=2))
  })
  mrres
}

p1 <- mr_scatter_plot(res, dat_clumped)
p1 <- p1[[1]] + pretty_plot(fontsize=6) + L_border() + theme(legend.position ="none") +
  geom_hline(yintercept=0,linetype="dashed",size=0.35)
p1

res_single <- mr_singlesnp(dat_clumped,single_method="mr_meta_fixed",
                           all_method=methods_to_use)
p2 <- mr_forest_plot(res_single)
p2 <- p2[[1]]+ pretty_plot(fontsize=6) + L_border() + theme(legend.position = "none")

res_loo <- mr_leaveoneout(dat_clumped,method = TwoSampleMR::mr_ivw)
p3 <- mr_leaveoneout_plot(res_loo)
p3 <- p3[[1]]+ pretty_plot(fontsize=6) + L_border() + theme(legend.position = "none")

# Save plots 
if (TRUE){
  cowplot::ggsave2(p1, file="../output/telomere_length/mendelian_randomization/MR_scatterplot.pdf", 
                   width=3.3, height=3.3,units="cm")
  cowplot::ggsave2(p2, file="../output/telomere_length/mendelian_randomization/MR_singleSNP.pdf", 
                   width=2.5, height=3)
  cowplot::ggsave2(p3, file="../output/telomere_length/mendelian_randomization/MR_leave_one_out.pdf", 
                   width=2.5, height=3)
}