library(data.table)
library(tidyverse)
library(BuenColors)

sentinels <- fread("../output/replication_joint_statistics/MVP.mpn_expanded.hasJAK2.tsv")

# Assume overall familial RR to first-degree relatives of cases
overall_RR = 4.93

# Formula taken from Michailidou et al., 2017 Nature
get_heritability <- function(sentinels){
  variant_rr <- sapply(seq(1,nrow(sentinels)),function(i){
    p <- sentinels[i,]$maf
    beta <- sentinels[i,]$effect
    se <- sentinels[i,]$stderr
    lambda <- overall_RR
    
    rr <- p *(1-p)*(beta^2-se^2)/log(lambda)
    
    return(rr)
    # ((1/log(overall_RR))*(1/rr)*(2*p*q*(r-1)/(p*r+q)^3))^2
  })
  return(sum(variant_rr))
}

get_heritability(sentinels)
get_heritability(sentinels %>% filter(meta_pval < 5e-8))