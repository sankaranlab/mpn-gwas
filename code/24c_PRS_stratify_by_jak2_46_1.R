library(tidyverse)
library(data.table)
library(pROC)
library(plotROC)
library(rms)
library(scales)
library(BuenColors)
library(forestplot)

# JAK2 46/1
carrier_path <- "JAK2 46/1 status"
jak2_461 <- fread(carrier_path)
colnames(jak2_461)[2] <- "jak2_461"
jak2_461$IID <- as.character(jak2_461$IID)
other_covariates <- merge(other_covariates,jak2_461,by="IID")

# Read in other covariates: age, sex, genotyping array, top PCs
covariates_path <- "filepath"
other_covariates <- fread(covariates_path)
other_covariates$IID <- as.character(other_covariates$IID)

# Read in PRS scores
toplot <- T
score_paths <- "file path for PRS scores"
scalar <- "1e-5"
scores <- fread(cmd=paste0(score_paths,scalar,".ALL.sscore.gz"))

ids <- scores$`#IID`
scores <- scores[,-1] %>% data.matrix() %>% rowSums()
scores.df <- data.frame(IID = ids, PRS= scores)
scores.df$IID <- as.character(scores.df$IID)

# Add PRS scores to the other covariates
all_covariates <- merge(other_covariates, scores.df,by="IID") %>% dplyr::rename(pheno="y")

# Calculate percentiles 
all_covariates$percentile_PRS <- cut(all_covariates$PRS, quantile(all_covariates$PRS, prob = 0:100 / 100, names = T),
                                     include.lowest=TRUE,labels=seq(1,100)) %>% as.integer()
all_covariates$decile_PRS <- cut(all_covariates$PRS, quantile(all_covariates$PRS, prob = 0:10 / 10, names = T),
                                 include.lowest=TRUE,labels=seq(1,10)) %>% as.integer()

# Compare risk stratified between carrier and PRS
all_covariates <- all_covariates %>% mutate(category = ifelse(jak2_461 > 0.1, "carrier","noncarrier"),
                                            PRS_category = case_when(decile_PRS <= 2 ~ "low",
                                                                     decile_PRS >= 8 ~ "high",
                                                                     TRUE ~ "intermediate"))
all_covariates <- all_covariates %>% mutate(PRS_carrier_category = paste(category,PRS_category,sep="-"))

# Divide into JAK2 46/1 carriers vs. non-carriers
carriers <- all_covariates %>% filter(category == "carrier")
noncarriers <- all_covariates %>% filter(category == "noncarrier")
dim(carriers)

table(all_covariates$pheno)
table(noncarriers$PRS_category,noncarriers$pheno)[,1]

# Forest plot
category_regression <- function(df,cat,baseline="noncarrier-intermediate"){
  print(cat)
  subset <- df %>% filter(PRS_carrier_category %in% c(baseline,cat)) %>% 
    mutate(PRS_carrier_category = ifelse(PRS_carrier_category == baseline, "baseline","risk"))
  covars_to_use <- subset %>% dplyr::select(age,sex,genotyping.array,paste0("PC",seq(1,10)),PRS_carrier_category)
  pheno <- subset$pheno
  
  prevalence <- 100*nrow(subset[subset$pheno == 1,]) / nrow(subset) 
  
  if (cat == baseline){
    d <- data.frame(category=cat,prevalence = prevalence,OR=1,bottom95CI=1,top95CI=1,p=NA)
  } else{
    results <- glm(pheno ~ ., data=covars_to_use,family="binomial") %>% summary()
    beta <- results$coefficients[nrow(results$coefficients),"Estimate"]
    error <- results$coefficients[nrow(results$coefficients),"Std. Error"]
    pval <- results$coefficients[nrow(results$coefficients),"Pr(>|z|)"]
    OR <- exp(beta)
    top95ci <- exp(beta + 1.96*error)
    bottom95ci <- exp(beta - 1.96*error)
    
    d <- data.frame(category=cat,
                    prevalence = prevalence,
                    OR=round(OR,3),
                    bottom95CI=round(bottom95ci,3),
                    top95CI=round(top95ci,3),
                    p=pval)
  }
  return(d)
}

ORs <- lapply(unique(all_covariates$PRS_carrier_category),
       function(y) category_regression(df=all_covariates,cat=y,baseline="noncarrier-intermediate")) %>% 
  bind_rows   %>% arrange(category)
ORs <- ORs %>% mutate(p = ifelse(p<0.0001,"<0.0001",round(p,3)))
ORs

tabletext<-cbind(
  c("JAK2 46/1", "Carrier", "Carrier", 
    "Carrier", "Noncarrier", "Noncarrier", "Noncarrier"),
  c("Polygenic Score", rep(c("High","Intermediate","Low"),2)),
  c("N of 1086 cases",table(carriers$PRS_category,carriers$pheno)[,2],table(noncarriers$PRS_category,noncarriers$pheno)[,2]),
  c("N of 407155 controls",table(carriers$PRS_category,carriers$pheno)[,1],table(noncarriers$PRS_category,noncarriers$pheno)[,1]))
tabletext<- cbind(tabletext,c("OR [95% CI]",paste0(ORs$OR," [", ORs$bottom95CI,"; ", ORs$top95CI,"]")))
tabletext<- cbind(tabletext,c("p-value",ORs$p))
tabletext


pdf(file="../output/PRS/PRS_JAK2_carrier_forestplot.pdf", width = 7.5, height = 2)  
par(cex.main=0.8,mar=c(1,1,1,1))

forestplot(tabletext,
           mean = c(NA,ORs$OR),
           lower = c(NA,ORs$bottom95CI),upper = c(NA,ORs$top95CI),
           graph.pos=5,
           zero=1,
           ci.vertices = TRUE, ci.vertices.height = 0.05,
           col=fpColors(box="black",line="black"),
           txt_gp = fpTxtGp(label=gpar(cex=0.5),
                            ticks=gpar(cex=0.5),
                            xlab=gpar(cex=0.5)),
           line.margin=unit(0,"mm"),
           boxsize=0.1)
dev.off()
