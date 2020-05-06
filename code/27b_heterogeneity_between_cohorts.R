library(data.table)
library(tidyverse)
library(BuenColors)
library(cowplot)
library(meta)
library(getmstatistic)  # for calculating M statistics
library(gridExtra)      # for generating tables

# Read in discovery sumstats
sentinels <- fread("../data/abf_finemap/MPN_arraycovar_meta_finngen_r4_abf_cojo_PP0.001_annotated.bed") %>% 
  filter(sentinel=="yes")
genomewide <- fread("../data/meta-gwas/sentinels/MPN_arraycovar_meta_finngen_r4_gcta_cojo_combined.5e-8.txt")
genomewide_vars <- sentinels %>% filter(rsid %in% genomewide$SNP)

ind_cohort_sumstats <- fread("individual_cohort_sumstats_file")
ind_cohort_sentinels <- ind_cohort_sumstats %>% filter(UKID %in% sentinels$var)

mvp_icd <- ind_cohort_sentinels %>% filter(cohort != "MVP_jak2")
mvp_jak2 <- ind_cohort_sentinels %>% filter(cohort != "MVP_jak2_or_mpn")
  
# Run M analysis on all sentinel variants for MVP ICD phenotype
mvp_icd_Mstatistic <- getmstatistic(mvp_icd$beta, 
                                       mvp_icd$se, 
                                       mvp_icd$UKID, 
                                       mvp_icd$cohort)
mvp_icd_Mstatistic

# Run M analysis on all sentinel variant for MVP JAK2 V617F phenotype
mvp_jak2_Mstatistic <- getmstatistic(mvp_jak2$beta, 
                                       mvp_jak2$se, 
                                       mvp_jak2$UKID, 
                                       mvp_jak2$cohort)
mvp_jak2_Mstatistic
dframe <- mvp_jak2_Mstatistic$M_dataset

# Run M analysis on all sentinel variant for all phenos
mvp_jak2_Mstatistic <- getmstatistic(ind_cohort_sentinels$beta, 
                                     ind_cohort_sentinels$se, 
                                     ind_cohort_sentinels$UKID, 
                                     ind_cohort_sentinels$cohort)
mvp_jak2_Mstatistic
dframe <- mvp_jak2_Mstatistic$M_dataset

# Write output table
dframe %>% filter(variant_names_in == "11:108143456_C_G")
output <- dframe %>% dplyr::select(study_names_in,variant_names_in, beta_in, lambda_se_in, 
                                   M, M_se, bonfpvalue, tau2, I2, Q) %>% unique()
fwrite(output,file="../output/replication_joint_statistics/between_cohort_heterogeneity_all.tsv",sep="\t")

# Write text file of sentinel UKIDs 
if (FALSE){
  sentinels <- fread("../data/abf_finemap/MPN_arraycovar_meta_finngen_r4_abf_cojo_PP0.001_annotated.bed") %>% filter(sentinel=="yes")
  fwrite(as.data.frame(sentinels$var),file="../sentinel_UKIDs.txt",sep="\t",col.names=FALSE)
}