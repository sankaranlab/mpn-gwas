library(tidyverse)
library(data.table)
library(pROC)
library(plotROC)
library(rms)
library(scales)
library(BuenColors)

# Read in other covariates: age, sex, genotyping array, top PCs
covariates_path <- "filepath"
other_covariates <- fread(covariates_path)
other_covariates$IID <- as.character(other_covariates$IID)

# Read in PRS scores
prune_and_thresh <- T
toplot <- T
thresholds <- c("1","5e-3","1e-3","5e-4","1e-4","5e-5","1e-5","1e-6","5e-8")

# Get AUROC for all p thresholds
all_scores <- lapply(thresholds,function(scalar){
  print(scalar)
  score_paths <- "file path for PRS scores"
  scores <- fread(paste0(score_paths,scalar,".ALL.sscore.gz"))
  scored_vars <- fread(paste0(score_paths,scalar,".ALL.scored.vars"),header=FALSE)
  clumped_vars <- fread(paste0(score_paths,scalar,".ALL.clumped"),header=FALSE)
  
  ids <- scores$`#IID`
  scores <- scores[,-1] %>% data.matrix() %>% rowSums()
  scores.df <- data.frame(IID = ids, PRS= scores,p_thresh=scalar)
  scores.df$IID <- as.character(scores.df$IID)
  
  # Add PRS scores to the other covariates
  all_covariates <- merge(other_covariates, scores.df,by="IID") %>% dplyr::rename(pheno="y")
  
  pheno <- all_covariates$pheno
  covars_to_use <- all_covariates %>% dplyr::select(age,sex,genotyping.array,paste0("PC",seq(1,10)),PRS)
  
  # Perform logistic regression
  prs_model <- glm(pheno ~ ., data=covars_to_use,family="binomial")
  predpr <- predict(prs_model,type=c("response"))
  
  # ROC
  results <- roc(pheno,predpr)
  
  output <- data.frame(p_thresh = scalar,
                       num_variants = paste(nrow(scored_vars),"/",nrow(clumped_vars)),
                       auroc = as.numeric(results$auc))
  return(output)
}) %>% bind_rows()

all_scores

if (TRUE){
  write.table(all_scores,file="../output/PRS/r4_PRS_AUROC_varying_p_thresholds.tsv",
              quote = FALSE, sep = "\t", col.names = TRUE, row.names = F)
}