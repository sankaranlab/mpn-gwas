library(tidyverse)
library(data.table)
library(pROC)
library(ROCR)
library(plotROC)
library(rms)
library(scales)
library(BuenColors)

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

# Restrict by a certain age
if (FALSE){
  all_covariates <- all_covariates %>% filter(age > 60)
}

# Calculate percentiles 
all_covariates$percentile_PRS <- cut(all_covariates$PRS, quantile(all_covariates$PRS, prob = 0:100 / 100, names = T),
                                     include.lowest=TRUE,labels=seq(1,100)) %>% as.integer()
all_covariates$decile_PRS <- cut(all_covariates$PRS, quantile(all_covariates$PRS, prob = 0:10 / 10, names = T),
                                 include.lowest=TRUE,labels=seq(1,10)) %>% as.integer()

pheno <- all_covariates$pheno
covars_to_use <- all_covariates %>% dplyr::select(age,sex,genotyping.array,paste0("PC",seq(1,10)),PRS)

# Perform logistic regression
prs_model <- glm(pheno ~ ., data=covars_to_use,family="binomial")
predpr <- predict(prs_model,type=c("response"))

# ROC
roc(pheno,predpr)

no_prs_model <- glm(pheno ~ ., data=covars_to_use %>% dplyr::select(-PRS),family="binomial")
noPRS_predpr <- predict(no_prs_model,type=c("response"))
roc(pheno,noPRS_predpr)

roc_df <- rbind(
  data.frame(pred = predpr,group = "PRS",y = pheno),
  data.frame(pred = noPRS_predpr,group = "no PRS", y = pheno)
)

# ROC curve
roc_plot <- ggplot(roc_df, aes(d = y, m = pred, color = group)) +
  geom_roc(n.cuts = 0) + pretty_plot(fontsize= 7) + L_border() +
  labs(x = "1 - Specificity", y = "Sensitivity", color = "") +
  scale_color_manual(values=jdb_palette("brewer_spectra")[c(7,1)]) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  theme(legend.position = "none") + 
  coord_fixed()

# Calculate Nagelkerke's pseudo-R2 to estimate explained variance
R2_difference <- function(variable,pheno,covars){
  diff <- lrm(pheno ~ .,data=covars)$stats["R2"] - lrm(pheno ~ .,data=covars %>% dplyr::select(-variable))$stats["R2"]
  df <- data.frame(variable = paste(variable,collapse="/"),incremental_r2 = diff)
  return(df)
}

lrm(pheno ~ .,data=covars_to_use)$stats["R2"]
lrm(pheno ~ .,data=covars_to_use %>% dplyr::select(-PRS))$stats["R2"]

variables <- list("PRS","sex","age","genotyping.array",paste0("PC",seq(10)))
incremental_r2 <- lapply(variables,function(y) R2_difference(variable = y, pheno =pheno,covars=covars_to_use)) %>% bind_rows()
incremental_r2$variable[5] <- "top 10 PCs"
incremental_r2$variable <- factor(incremental_r2$variable, levels=unique(incremental_r2$variable))

env_covars <- all_covariates %>% dplyr::select(age,sex,genotyping.array,paste0("PC",seq(1,10)),irradiation)
variables <- list("irradiation")
env_incremental_r2 <- lapply(variables,function(y) R2_difference(variable = y, pheno =pheno,covars=env_covars)) %>% bind_rows()


p1 <- ggplot(incremental_r2,aes(x=variable,y=100*incremental_r2)) +
  geom_bar(position=position_dodge(), stat="identity", fill = "firebrick", color = "black")  + 
  labs(x="",y="Incremental pseudo-R2 (%)")+
  pretty_plot(fontsize = 7) + L_border()+
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust=1)) 
p1 

# Calculate ORs for different deciles of PRS 
# Compare individuals from decile X to those from baseline decile
decile_regression <- function(decile,baseline=1){
  #stopifnot(decile != 1)
  print(decile)
  
  subset <- all_covariates %>% filter(decile_PRS %in% c(baseline,decile)) %>% mutate(decile_PRS = ifelse(decile_PRS == baseline, "baseline","risk"))
  covars_to_use <- subset %>% dplyr::select(age,sex,genotyping.array,paste0("PC",seq(1,10)),decile_PRS)
  pheno <- subset$pheno
  
  prevalence <- 100*nrow(subset[subset$pheno == 1,]) / nrow(subset) 
  
  if (decile ==baseline){
    d <- data.frame(decile=decile,n = nrow(subset[subset$pheno == 1,]),prevalence = prevalence,OR=1,bottom95CI=1,top95CI=1)
  } else{
    results <- glm(pheno ~ ., data=covars_to_use,family="binomial") %>% summary()
    beta <- results$coefficients[nrow(results$coefficients),"Estimate"]
    error <- results$coefficients[nrow(results$coefficients),"Std. Error"]
    OR <- exp(beta)
    top95ci <- exp(beta + 1.96*error)
    bottom95ci <- exp(beta - 1.96*error)
    
    d <- data.frame(decile=decile,
                    n = nrow(subset[subset$pheno == 1,]),
                    prevalence = prevalence,
                    OR=round(OR,3),
                    bottom95CI=round(bottom95ci,3),
                    top95CI=round(top95ci,3))
  }
  return(d)
}

decile_ORs <- lapply(seq(1,10),function(y) decile_regression(y,baseline=5)) %>% bind_rows 
decile_ORs
decile_ORs_compared_to_1 <- lapply(seq(1,10),function(y) decile_regression(y,baseline=1)) %>% bind_rows 
decile_ORs_compared_to_1

p2 <- ggplot(decile_ORs,aes(x=decile,y=OR)) +
  geom_point(shape = 15) + 
  geom_errorbar(aes(ymin = bottom95CI, ymax = top95CI), colour="black", width=.1) +
  geom_hline(yintercept = 1, linetype = 2) +
  labs(x="Decile of PRS (1 = lowest)",y="Odds Ratio") +
  scale_x_continuous(breaks=seq(1,10)) +
  scale_y_continuous(breaks=seq(1,6)) +
  theme(legend.position="none") +
  pretty_plot(fontsize = 7) + L_border() 
p2

p3 <- ggplot(decile_ORs_compared_to_1,aes(x=decile,y=OR)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(ymin = bottom95CI, ymax = top95CI), colour="black", width=.1) +
  geom_hline(yintercept = 1, linetype = 2) +
  labs(x="Decile of PRS (1 = lowest (reference))",y="Odds Ratio") +
  scale_x_continuous(breaks=seq(1,10)) +
  scale_y_continuous(breaks=seq(1,6)) +
  theme(legend.position="none") +
  pretty_plot(fontsize = 7) + L_border() 

p4 <- ggplot(decile_ORs,aes(x=decile,y=prevalence)) +
  geom_point(shape = 15) + 
  labs(x="Decile of PRS (1 = lowest (reference))",y="Prevalence of MPN (%)") +
  scale_x_continuous(breaks=seq(1,10)) +
  theme(legend.position="none") +
  pretty_plot(fontsize = 7) + L_border() 

if (toplot){
  cowplot::ggsave2(roc_plot,file=paste0("../output/PRS/r4_PRS.pruned.",scalar,"_roc_comparison.pdf"),width=1.5,height=1.5)  
  cowplot::ggsave2(p1,file=paste0("../output/PRS/r4_PRS.pruned.",scalar,"_incremental_r2_comparison.pdf"),width=1.5,height=1.5)  
  cowplot::ggsave2(p2, file=paste0("../output/PRS/r4_PRS.pruned.",scalar,"_decile_ORs.pdf"),width = 1.5, height=1.5)
  cowplot::ggsave2(p3, file=paste0("../output/PRS/r4_PRS.pruned.",scalar,"_decile_ORs_compared_to_1.pdf"),width = 1.5, height=1.5)
  cowplot::ggsave2(p4, file=paste0("../output/PRS/r4_PRS.pruned.",scalar,"_decile_prevalence.pdf"),width = 1.5, height=1.5)
}


if (FALSE){
  write.table(all_covariates[,-1],file=paste0("../output/PRS/r4_PRS_scores_covariates_PT_",scalar,".txt"),quote = FALSE, sep = "\t", col.names = TRUE, row.names = F)
}