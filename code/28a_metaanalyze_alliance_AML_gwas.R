library(data.table)
library(tidyverse)

# Read in Alliance sumstats
alliance <- fread("../data/other_myeloid_neoplasms/alliance/MPN_risk_loci_p1e5inAllianceAML.txt")

# Verify that alliance sub-cohorts are oriented
alliance %>% filter(USA1.EA != USA2.EA)

# MPN stats
CS.df <- fread("../data/abf_finemap/MPN_arraycovar_meta_finngen_r4_abf_cojo_PP0.001_annotated.bed") 

# Inverse variance weighted meta-analysis function
ivw_meta <- function(betas,ses){
  weights <- 1/(ses^2)
  
  m_se <- sqrt(1/(sum(weights)))
  m_beta <- sum(betas*weights) / sum(weights)
  m_z <- m_beta / m_se
  
  m_pval <- 2*pnorm(-abs(m_z))
  return(c(meta_beta = m_beta, meta_se=  m_se,meta_pval = m_pval))
}

combined_cases <- 1183
lapply(alliance$UKID,function(var){
  subset <- alliance %>% dplyr::filter(UKID == var)
  
  # Construct ID based on chr:pos_OA_EA format
  test_UKID <- alliance %>% dplyr::filter(UKID == var) %>%
    mutate(UKID_2 = paste0(CHR,":",POS,"_",USA1.OA,"_",USA1.EA)) %>% pull(UKID_2)
  
  meta_stats <- ivw_meta(betas = c(subset$USA1.Beta,subset$USA2.Beta),ses = c(subset$USA1.SE,subset$USA2.SE))
  
  out <- data.frame(trait = "AML",
                    UKID = var,
                    UKID_2 = test_UKID,
                    n_cases = combined_cases,
                    beta = meta_stats["meta_beta"],
                    se = meta_stats["meta_se"],
                    pval = meta_stats["meta_pval"])
  out
}) %>% bind_rows -> out

# Orient with respect to UKBB
out %>% filter(UKID != UKID_2)
oriented <- out %>% 
  mutate(UKID_true = ifelse(UKID %in% CS.df$var,UKID,UKID_2),
         beta_true = ifelse(UKID %in% CS.df$var,beta,-1*beta)) %>% 
  dplyr::select(-c(UKID,UKID_2,beta)) %>%
  dplyr::rename(UKID = "UKID_true",beta= "beta_true") %>% dplyr::select(UKID, n_cases,beta,se,pval)

# Write meta-analyzed sumstats
fwrite(out, file="../data/other_myeloid_neoplasms/alliance/MPN_risk_loci_Alliance_meta.txt",sep="\t")
