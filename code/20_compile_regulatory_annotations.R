library(tidyverse)
library(data.table)
library(BuenColors)
library(GenomicRanges)
library(readxl)
library(rtracklayer)
library(GenomicRanges)
library(magrittr)
"%ni%" <- Negate("%in%")

trait <- "MPN_arraycovar_meta_finngen_r4"
# Read in CS
CS.df <- fread(paste0("../data/abf_finemap/",trait,"_abf_cojo_PP0.001_annotated.bed") )%>% 
  mutate(ref = str_split_fixed(var,"_",3)[,2],alt = str_split_fixed(var,"_",3)[,3]) %>%
  group_by(var) %>% dplyr::slice(which.max(PP)) %>% mutate(pos = start) %>% ungroup()

# Read in VEP
CS.df.vep <- fread(paste0("../output/VEP/",trait,".CS_PP001_VEP_annotated.tsv"))

# Merge with regulomeDB
regulomeDB <- fread("../output/regulatory_annotations/MPN_r4_PP001_regulomedb_scores_only.txt")
CS.df.vep <- merge(CS.df.vep,regulomeDB,by="rsid")
CS.df.vep.PP01 <- CS.df.vep %>% filter(PP>0.01)

# Merge with hemeATAC
peaksdf <- fread("../data/atac/26August2017_EJCsamples_allReads_250bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
idx <- findOverlaps(peaks,GRanges(CS.df.vep))
CS.df.vep <- CS.df.vep %>% mutate(hemeATAC = ifelse(var %in% CS.df.vep[idx@to,]$var,"yes","no"))

# Add gene targets
variant_to_gene <- fread("../output/target_genes/r4_scored_noABC_target_genes_variant_list_0.01.tsv") %>%
  dplyr::select(var,gene) %>% group_by(var) %>% 
  summarize(target_genes =paste(gene, collapse = ','))
CS.df.vep <- left_join(CS.df.vep,variant_to_gene[,c("var","target_genes")],by="var")

# Add TF motifs
motif_creations <- fread("../output/transcription_factors/r4_mpn_motifbreakr_PPmerged.txt") %>% 
  filter(alleleAlt > alleleRef) %>%  dplyr::select(var,geneSymbol) %>% group_by(var) %>%
  summarize(motif_creations =paste(geneSymbol, collapse = ','))
motif_disruptions <- fread("../output/transcription_factors/r4_mpn_motifbreakr_PPmerged.txt") %>% 
  filter(alleleAlt < alleleRef) %>%  dplyr::select(var,geneSymbol) %>% group_by(var) %>%
  summarize(motif_disruptions =paste(geneSymbol, collapse = ','))
CS.df.vep <- left_join(CS.df.vep,motif_creations[,c("var","motif_creations")],by="var")
CS.df.vep <- left_join(CS.df.vep,motif_disruptions[,c("var","motif_disruptions")],by="var")

CS.df.vep[is.na(CS.df.vep)] <- ""

# Import PhyloP
phyloP_file <- "/Volumes/broad_sankaranlab/tools/conservation/hg19.100way.phyloP100way.bw"
pp <- rtracklayer::import.bw(phyloP_file, which = GRanges(CS.df.vep), as = "NumericList")
CS.df.vep$Phylop_100way <- sapply(pp, function(x) x[1])

# Output supplementary table
annotated <- CS.df.vep %>% filter(PP>0.001) %>% 
  dplyr::select(rsid, var, region, sentinel, seqnames,PP,effect,stderr,pvalue,maf,target_genes,
                Phylop_100way,CADD_PHRED,most_severe_consequence,rdb,hemeATAC,motif_disruptions,motif_creations) %>%
  arrange(region,desc(PP))
fwrite(annotated,file="../output/regulatory_annotations/finemap_vars_with_annotations_PP001.tsv",sep="\t")












