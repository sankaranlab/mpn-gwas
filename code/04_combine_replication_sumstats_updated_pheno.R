library(data.table)
library(tidyverse)
library(BuenColors)
library(cowplot)

# Read in discovery sumstats
CS.df <- fread("../data/abf_finemap/MPN_arraycovar_meta_finngen_r4_abf_cojo_PP0.001_annotated.bed")
sentinels <- CS.df %>% filter(sentinel=="yes")

# Read in replication statistics and reorient
name <- "mpn_expanded.hasJAK2"

if (FALSE){
  replication <- fread(paste0("../data/MVP_replication/revised_MVP/",name,".glm.logistic.txt")) %>% unique()

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

  # Take minimum P value for duplicate variant rows
  replication <- replication %>% group_by(UKID) %>% dplyr::slice(which.min(P)) %>% ungroup()
  
  fwrite(replication,file=paste0("../data/MVP_replication/revised_MVP/",name,".glm.logistic.oriented.txt"),sep="\t")
}

replication <- fread(paste0("../data/MVP_replication/revised_MVP/",name,".glm.logistic.oriented.txt"))
merged <- merge(CS.df,replication,by.x="var",by.y="UKID")

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

merged <- merged %>% group_by(region) %>% dplyr::slice(which.min(meta_pval)) %>% ungroup() %>% as.data.frame()
merged %>% dplyr::select(var,rsid,pvalue,P,meta_pval,region_rank) %>% arrange(meta_pval)

# Binomial test for directional concordance
merged <- merged %>% mutate(direction = ifelse((effect > 1 & beta >1 )| (effect < 1 & beta <1),1,-1))
merged %>% dplyr::select(effect,beta,direction)
directions <- as.numeric(table(merged$direction))
if (length(directions)==1){directions[2] <- directions[1]; directions[1]<-0}

binom.test(directions[2], nrow(merged), p = 0.5, alternative = "greater", conf.level = 0.95)

# How many significant at 5% level
nominal <- merged %>% filter(P < 0.05)
binom.test(nrow(nominal), nrow(merged), p = 0.05, alternative = "greater", conf.level = 0.95)

# Plot effect size correlation
summary(lm(effect~beta,merged))

p1 <- ggplot(merged,aes(x=effect,y=beta))+
  geom_point()+
  geom_smooth(method = "lm",formula=y~x)+
  labs(x="Discovery log(OR)",y="Replication log(OR)") +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept=0, linetype="dashed")+
  pretty_plot(fontsize = 7) + L_border()
p1
cowplot::ggsave2(p1,filename="../output/replication_joint_statistics/mpn_replication_logORs.pdf",
                 height=2.5,width=2.5)
cor.test(merged$effect,merged$beta,method = "pearson")

merged <- merged %>% arrange(pvalue)
merged %>% filter(meta_pval < 5e-8,pvalue > 5e-8) 
merged %>% filter(meta_pval > 5e-8,pvalue < 5e-8) 


merged %>% dplyr::select(var,pvalue,P,fishers_p) %>% filter(fishers_p < 5e-8,pvalue > 5e-8) 

write.table(merged,file=paste0("../output/replication_joint_statistics/MVP.",name,".tsv"),
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = F)

