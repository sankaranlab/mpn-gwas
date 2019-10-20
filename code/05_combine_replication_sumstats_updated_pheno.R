library(data.table)
library(tidyverse)
library(BuenColors)
library(cowplot)

# Read in discovery sumstats
CS.df <- fread("../data/abf_finemap/MPN_CML_abf_cojo_PP0.001_annotated.bed") %>% filter(sentinel=="yes")

# Read in replication statistics and reorient and save
name <- "mpnRepLoci.hasJAK2"
if (FALSE){
  replication <- fread(paste0("../data/MVP_replication/",name,".glm.logistic.txt"))

  replication <- replication %>% drop_na() %>%
    dplyr::rename(CHR = "#CHROM") %>% 
    mutate(UKID = paste0(CHR,":",POS,"_",REF,"_",ALT),UKID_2 = paste0(CHR,":",POS,"_",ALT,"_",REF)) %>%
    mutate(beta = log(OR), se = abs(log(OR)/T_STAT)) 
  # beta is log(OR) and the SE is abs(log(OR)/STAT)
  
  # Orient UKID in the same way
  replication <- replication %>% 
    mutate(UKID_true = ifelse(UKID %in% CS.df$var,UKID,UKID_2),
           beta_true = ifelse(UKID %in% CS.df$var,beta,-1*beta),
           ALT_true = ifelse(UKID %in% CS.df$var,ALT,REF),
           REF_true = ifelse(UKID %in% CS.df$var,REF,ALT)) %>% 
    dplyr::select(-c(UKID,UKID_2,beta,REF,ALT)) %>%
    dplyr::rename(UKID = "UKID_true",beta= "beta_true",REF = "REF_true",ALT="ALT_true") %>% 
    mutate(OR = exp(beta)) 
  fwrite(replication,file=paste0("../data/MVP_replication/",name,".glm.logistic.oriented.txt"),sep="\t")
}

replication <- fread(paste0("../data/MVP_replication/",name,".glm.logistic.oriented.txt"))
merged <- merge(CS.df,replication,by.x="var",by.y="UKID")

# Binomial test for directional concordance
merged <- merged %>% mutate(direction = ifelse((effect > 1 & beta >1 )| (effect < 1 & beta <1),1,-1))
merged %>% dplyr::select(effect,beta,direction)
directions <- as.numeric(table(merged$direction))

binom.test(directions[2], nrow(merged), p = 0.5, alternative = "greater", conf.level = 0.95)

# How many significant at 5% level
nominal <- merged %>% filter(P < 0.05)
binom.test(nrow(nominal), nrow(merged), p = 0.05, alternative = "greater", conf.level = 0.95)

# Plot effect size correlation
summary(lm(effect~beta,merged))

ggplot(merged,aes(x=beta,y=effect))+
  geom_point()+
  pretty_plot(fontsize = 8) + L_border()

# Inverse variance weighted meta-analysis
ivw_meta <- function(betas,ses){
  weights <- 1/(ses^2)
  
  meta_se <- sqrt(1/(sum(weights)))
  meta_beta <- sum(betas*weights) / sum(weights)
  meta_z <- meta_beta / meta_se
  
  meta_pval <- 2*pnorm(-abs(meta_z))
  return(meta_pval)
}
merged$meta_pval <- sapply(merged$var,function(UKID){
  subset <- merged %>% filter(var == UKID)
  ivw_meta(betas = c(subset$effect,subset$beta),ses = c(subset$stderr,subset$se))
})

merged %>% dplyr::select(var,pvalue,P,meta_pval) %>% arrange(meta_pval)

write.table(merged,file="../output/replication_joint_statistics/MVP_replication_mvp_JAK2_only.tsv",
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = F)

