library(data.table)
library(tidyverse)
library(BuenColors)
library(cowplot)

# Gene annots
gencode_gr <- readRDS("../../bcx-finemap/data/annotations/gencode_filtered.rds")

# Read in sumstats
jak2 <- fread("../data/MVP_replication/revised_MVP/mpn_expanded.hasJAK2.glm.logistic.oriented.txt")
jak2_or_mpn <- bind_rows(fread("../data/MVP_replication/mpnRepLoci-phe4.hasMPNorJAK2.glm.logistic.oriented.txt"),
                         fread("../data/MVP_replication/mpnRep-chr13Loci-phe4.hasMPNorJAK2.glm.logistic.oriented.txt"))
setdiff(jak2$UKID,jak2_or_mpn$UKID)

CS.df <- fread("../data/abf_finemap/MPN_arraycovar_meta_finngen_r4_abf_cojo_PP0.001_annotated.bed") 
for_extraction <- CS.df %>% filter(PP>0.05) 
sentinels <- fread("../output/significant_loci_tables/MPN_r4_genomewide_loci.tsv") %>% dplyr::rename(rsid = "SNP")
genomewide <- fread("../data/meta-gwas/sentinels/MPN_arraycovar_meta_finngen_r4_gcta_cojo_combined.5e-8.txt")
genomewide_vars <- sentinels %>% filter(rsid %in% genomewide$SNP)

# Pick select rsids to display
select_vars <- CS.df %>% filter(rsid %in% c("rs7868130","rs7705526","rs2853677"))
common_vars <- intersect(jak2$UKID,sentinels$var)

sumstats_file <- "file path for cohort-specific sum stats"
cohort_sumstats <- fread(sumstats_file) %>%
  dplyr::rename("beta" = BETA,"se" = SE, "P"=p.value,"REF"="Allele1","ALT"="Allele2")

tocombine <- jak2 %>% filter(UKID %in% cohort_sumstats$UKID) %>% mutate(cohort = "MVP_jak2") %>% 
  dplyr::select(cohort,UKID,REF,ALT,beta,se,P)
tocombine2 <- jak2_or_mpn %>% filter(UKID %in% cohort_sumstats$UKID) %>% mutate(cohort = "MVP_jak2_or_mpn") %>% 
  dplyr::select(cohort,UKID,REF,ALT,beta,se,P) %>% unique()

combined <- bind_rows(tocombine,
                      tocombine2,
                      cohort_sumstats %>% dplyr::select(cohort,UKID,REF,ALT,beta,se,P)) %>%
  arrange(P)
# Check that ref/alt alleles are oriented
combined %>% group_by(UKID) %>% summarise(distinct_alleles = n_distinct(REF)) %>% as.data.frame()

# Make UKBB direction positive
toflip <- combined %>% filter(cohort == "UKBB", beta < 0) %>% .$UKID
combined <- combined %>% mutate(beta=ifelse(UKID %in% toflip,-1*beta,beta))

# Add nearest gene
CS.df$nearest_gene <- gencode_gr[nearest(GRanges(CS.df),gencode_gr,ignore.strand=TRUE),]$gene_name

# Merge with RSIDs
combined <- merge(combined, CS.df[,c("var","rsid","nearest_gene","region")],by.x="UKID",by.y="var") %>%
  mutate(ID = paste0(rsid," (",nearest_gene,")")) %>% arrange(region)

# Convert BETA to ORs
combined$OR <- round(exp(combined$beta),2)
combined$OR_lowerCI <- round(exp(combined$beta - 1.96*combined$se),2)
combined$OR_upperCI <- round(exp(combined$beta + 1.96*combined$se),2)

# Factorize variables
combined$cohort <- factor(combined$cohort,levels=rev(c("UKBB","23andMe","FinnGen",
                                                   "MVP_jak2","MVP_jak2_or_mpn")))
combined$UKID <- factor(combined$UKID,levels=unique(combined$UKID))
combined$rsid <- factor(combined$rsid,levels=unique(combined$rsid))
combined$ID <- factor(combined$ID,levels=unique(combined$ID))
combined$sig <- ifelse(combined$P < 0.05,"yes","no")

# Write table of cohort sumstats
if (FALSE){
  fwrite(combined,file="../output/replication_joint_statistics/ind_cohort_sumstats.tsv",sep="\t")
}

# Graph
if (FALSE){
  selected <- ggplot(data=combined %>% filter(rsid %in% select_vars$rsid),
                     aes(x = cohort,y = OR,ymin = OR_lowerCI, ymax = OR_upperCI))+
    geom_pointrange(shape=15,color="black",fill="black",cex=0.25)+
    geom_hline(yintercept =1, linetype=2)+
    ylab("Odds Ratio")+ xlab("")+
    geom_errorbar(aes(ymin=OR_lowerCI, ymax=OR_upperCI),col="black",width=0.25,cex=0.5)+ 
    facet_wrap(~ID,strip.position="left",nrow=9,scales ="fixed") +
    coord_flip() + 
    pretty_plot(fontsize=6) 
  selected
  
  cowplot::ggsave2(selected,file="../output/locus_plots/ind_cohort_forestplot_selected.pdf",width=2,height=3)
  
  genomewide <- ggplot(data=combined %>% filter(rsid %in% sentinels$rsid),
                       aes(x = cohort,y = OR,ymin = OR_lowerCI, ymax = OR_upperCI))+
    geom_pointrange(shape=15,color="black",fill="black",cex=0.25)+
    geom_hline(yintercept =1, linetype=2)+
    ylab("Odds Ratio")+ xlab("")+
    geom_errorbar(aes(ymin=OR_lowerCI, ymax=OR_upperCI),col="black",width=0.25,cex=0.5)+ 
    facet_wrap(~ID,strip.position="left",ncol=3,scales ="free_x") +
    coord_flip() + 
    pretty_plot(fontsize=7)
  genomewide
  
  # Exclude JAK2_or_mpn
  combined_justJAK2 <- combined %>% filter(rsid %in% sentinels$rsid, cohort != "MVP_jak2_or_mpn") 
  combined_justJAK2$cohort <- as.character(combined_justJAK2$cohort)
  combined_justJAK2 <- combined_justJAK2 %>% mutate(cohort = ifelse(cohort == "MVP_jak2","MVP",cohort))
  combined_justJAK2$cohort <- factor(combined_justJAK2$cohort,levels=rev(c("UKBB","23andMe","FinnGen","MVP")))
  
  genomewide_justJAK2 <- ggplot(data=combined_justJAK2,
                       aes(x = cohort,y = OR,ymin = OR_lowerCI, ymax = OR_upperCI))+
    geom_pointrange(shape=15,color="black",fill="black",cex=0.25)+
    geom_hline(yintercept =1, linetype=2)+
    ylab("Odds Ratio")+ xlab("")+
    geom_errorbar(aes(ymin=OR_lowerCI, ymax=OR_upperCI),col="black",width=0.25,cex=0.5)+ 
    facet_wrap(~ID,strip.position="left",ncol=3,scales ="free_x") +
    coord_flip() + 
    pretty_plot(fontsize=7)
  genomewide_justJAK2
}

# Save
cowplot::ggsave2(genomewide,file="../output/locus_plots/ind_cohort_forestplot_genomewide.pdf",width=5,height=8)
cowplot::ggsave2(genomewide_justJAK2,file="../output/locus_plots/ind_cohort_forestplot_genomewide_justJAK2.pdf",width=5,height=8)
cowplot::ggsave2(p,file="../output/locus_plots/ind_cohort_forestplot_suggestive.pdf",width=6,height=8)
