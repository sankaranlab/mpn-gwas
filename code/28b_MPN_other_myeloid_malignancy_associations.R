library(tidyverse)
library(data.table)
library(BuenColors)

# Gene annots
gencode_gr <- readRDS("../data/annotations/gencode_gene_annotations.rds")

# Read in discovery sumstats
CS.df <- fread("../data/abf_finemap/MPN_arraycovar_meta_finngen_r4_abf_cojo_PP0.001_annotated.bed") 
# Add nearest gene
CS.df$nearest_gene <- gencode_gr[nearest(GRanges(CS.df),gencode_gr,ignore.strand=TRUE),]$gene_name
meta_genomewide <- fread("../output/significant_loci_tables/MPN_r4_genomewide_loci.tsv") %>% dplyr::rename(rsid = "SNP")
sentinels <- CS.df %>% filter(rsid %in% meta_genomewide$rsid)

# Traits and sample sizes
traits <- c("AML","MDS")
pheno_sizes <- fread("../data/other_myeloid_neoplasms/additional_pheno_casecontrol_numbers.txt")

# Inverse variance weighted meta-analysis function
ivw_meta <- function(betas,ses){
  weights <- 1/(ses^2)
  
  m_se <- sqrt(1/(sum(weights)))
  m_beta <- sum(betas*weights) / sum(weights)
  m_z <- m_beta / m_se
  
  m_pval <- 2*pnorm(-abs(m_z))
  return(c(meta_beta = m_beta, meta_se=  m_se,meta_pval = m_pval))
}

# Read in subtype sumstats and meta-analyze
lapply(traits, function(trait) {
  ukbb <- fread(paste0("../data/other_myeloid_neoplasms/ukbb/MPN_r4_sentinels_",trait,"_sumstats_ukbb.txt")) %>%
    dplyr::select(trait,UKID,SNPID,CHR,POS,Allele1,Allele2,BETA,SE,p.value) %>% 
    mutate(UKID_2 = paste0(CHR,":",POS,"_",Allele2,"_",Allele1))
  finngen <- fread(paste0("../data/other_myeloid_neoplasms/finngen/MPN_r4_sentinels_",trait,"_sumstats_finngen.txt")) %>%
    dplyr::select(trait,UKID,SNPID,chr,pos_hg19,ref,alt,beta,sebeta,pval) %>%
    mutate(UKID_2 = paste0(chr,":",pos_hg19,"_",alt,"_",ref))
  colnames(ukbb) <- colnames(finngen) <- c("trait","UKID","rsid","chr","pos","ref","alt","beta","sebeta","pval","UKID_2")
  
  # Orient stats in the same way
  ukbb_oriented <- ukbb %>% 
    mutate(UKID_true = ifelse(UKID %in% CS.df$var,UKID,UKID_2),
           beta_true = ifelse(UKID %in% CS.df$var,beta,-1*beta),
           alt_true = ifelse(UKID %in% CS.df$var,alt,ref),
           ref_true = ifelse(UKID %in% CS.df$var,ref,alt)) %>% 
    dplyr::select(-c(UKID,UKID_2,beta,ref,alt)) %>%
    dplyr::rename(UKID = "UKID_true",beta= "beta_true",REF = "ref_true",ALT="alt_true") %>% 
    mutate(OR = exp(beta)) 
  
  finngen_oriented <- finngen %>% 
    mutate(UKID_true = ifelse(UKID %in% CS.df$var,UKID,UKID_2),
           beta_true = ifelse(UKID %in% CS.df$var,beta,-1*beta),
           alt_true = ifelse(UKID %in% CS.df$var,alt,ref),
           ref_true = ifelse(UKID %in% CS.df$var,ref,alt)) %>% 
    dplyr::select(-c(UKID,UKID_2,beta,ref,alt)) %>%
    dplyr::rename(UKID = "UKID_true",beta= "beta_true",REF = "ref_true",ALT="alt_true") %>% 
    mutate(OR = exp(beta)) 
  
  # Get meta-pval
  common_vars <- intersect(ukbb_oriented$UKID,finngen_oriented$UKID)
  out_df <- lapply(unique(ukbb_oriented$UKID),function(var){
    # If variant in both UKBB and FinnGen, do IVW meta-analysis, otherwise just take UKBB stats
    if (var %in% common_vars){
      subset_1 <- ukbb_oriented %>% filter(UKID == var)
      subset_2 <- finngen_oriented %>% filter(UKID == var)
      
      meta_stats <- ivw_meta(betas = c(subset_1$effect,subset_2$beta),ses = c(subset_1$stderr,subset_2$se))
      
      cases  = sum(pheno_sizes %>% filter(Pheno == trait) %>% pull(Cases))
      controls  = sum(pheno_sizes %>% filter(Pheno == trait) %>% pull(Controls))
      
      out <- data.frame(trait = trait,
                        UKID = var,
                        n_cases = cases,
                        #n_controls = controls,
                        beta = meta_stats["meta_beta"],
                        se = meta_stats["meta_se"],
                        pval = meta_stats["meta_pval"])
      out
    } else{
      subset_1 <- ukbb_oriented %>% filter(UKID == var)
      
      cases  = pheno_sizes %>% filter(Pheno == trait,Cohort == "ukbb") %>% pull(Cases)
      controls  = pheno_sizes %>% filter(Pheno == trait,Cohort == "ukbb") %>% pull(Controls)
      
      out <- data.frame(trait = trait,
                        UKID = var,
                        n_cases = cases,
                        #n_controls = controls,            
                        beta = subset_1$beta,
                        se = subset_1$sebeta,
                        pval = subset_1$pval)
      out
    }
    
  }) %>% bind_rows()
  
  
  return(out_df)
}) %>% bind_rows() -> meta_myeloid_phenos

# Add on AML from Walker et al. ---------------------------------------------------------------
alliance_meta <- fread("../data/other_myeloid_neoplasms/alliance/MPN_risk_loci_Alliance_meta.txt")

common_vars <- intersect(alliance_meta$UKID,meta_myeloid_phenos$UKID)
lapply(common_vars,function(var){
  subset_1 <- alliance_meta %>% dplyr::filter(UKID == var)
  subset_2 <- meta_myeloid_phenos %>% dplyr::filter(UKID == var,trait == "AML")
  
  meta_stats <- ivw_meta(betas = c(subset_1$beta,subset_2$beta),ses = c(subset_1$se,subset_2$se))
  
  cases  = subset_1$n_cases + subset_2$n_cases
  
  out <- data.frame(trait = "AML",
                    UKID = var,
                    n_cases = cases,
                    beta = meta_stats["meta_beta"],
                    se = meta_stats["meta_se"],
                    pval = meta_stats["meta_pval"])
  out
}) %>% bind_rows -> out

# Other myeloid neoplasms
other_myeloid <- bind_rows(out, meta_myeloid_phenos %>% filter(trait == "MDS"))
other_myeloid %>% filter(pval < 0.05)

# Merge with MPN sumstats 
mpn <- sentinels %>% mutate(trait = "MPN",n_cases = 2949, n_controls = 835554) %>%
  dplyr::rename(UKID = "var",beta = "effect",se = "stderr",pval = "pvalue") %>%
  dplyr::select(trait,UKID,n_cases, beta,se,pval)

combined <- bind_rows(mpn,other_myeloid) 

# Merge with RSIDs
combined <- merge(combined, CS.df[,c("var","rsid","nearest_gene","region","maf")],by.x="UKID",by.y="var") %>%
  mutate(ID = paste0(rsid," (",nearest_gene,")")) %>% 
  arrange(region,factor(combined$trait,levels=c("MPN","AML","MDS")), desc(trait))
combined <- merge(combined, meta_genomewide[,c("rsid","risk")])

# Orient with respect to risk allele
combined <- combined %>% 
  mutate(beta_true = ifelse(str_split_fixed(combined$UKID,"_",3)[,3] == combined$risk,beta,-1*beta)) %>% 
  dplyr::select(-beta) %>%
  dplyr::rename(beta= "beta_true") 

# Convert BETA to ORs
combined$OR <- round(exp(combined$beta),2)
combined$OR_lowerCI <- round(exp(combined$beta - 1.96*combined$se),2)
combined$OR_upperCI <- round(exp(combined$beta + 1.96*combined$se),2)
combined$OR_95CI <- paste0(round(exp(combined$beta - 1.96*combined$se),2),"-",
                                  round(exp(combined$beta + 1.96*combined$se),2))

# Factorize variables
# combined$trait <- factor(combined$trait,levels=c("MPN","AML","MDS"))
combined$UKID <- factor(combined$UKID,levels=unique(combined$UKID))
combined$rsid <- factor(combined$rsid,levels=unique(combined$rsid))
combined$ID <- factor(combined$ID,levels=unique(combined$ID))
combined$sig <- ifelse(combined$pval < 0.05,"yes","no")

# Table
genomewide <- combined %>% filter(rsid %in% sentinels$rsid) %>% mutate()
subtype_comparison <- pivot_wider(genomewide,id_cols = c(rsid,nearest_gene,region,maf,risk),
                                  names_from=trait,values_from=c(pval,OR,OR_95CI,n_cases)) %>% 
  as.data.frame() %>% arrange(region) %>%
  dplyr::select(rsid,nearest_gene,risk,maf,pval_MPN,OR_MPN,OR_95CI_MPN,pval_AML,OR_AML,OR_95CI_AML)
head(subtype_comparison)

fwrite(subtype_comparison,file="../output/other_myeloid_neoplasms/other_myeloid_neoplasms_loci_comparison.tsv",sep="\t")

# Directional concordance
subtype_comparison <- fread("../output/other_myeloid_neoplasms/other_myeloid_neoplasms_loci_comparison.tsv") %>%
  mutate(direction = case_when(OR_MPN > 1 & OR_AML >1 ~ 1,
                               OR_MPN < 1 & OR_AML <1 ~ 1,
                               TRUE ~ -1))
directions <- as.numeric(table(subtype_comparison$direction))
if (length(directions)==1){directions[2] <- directions[1]; directions[1]<-0}

binom.test(directions[2], nrow(subtype_comparison), p = 0.5, alternative = "greater", conf.level = 0.95)

# How many significant at 5% level
nominal <- subtype_comparison %>% filter(direction == 1, pval_AML < 0.05)
binom.test(nrow(nominal), nrow(subtype_comparison), p = 0.05, alternative = "greater", conf.level = 0.95)


