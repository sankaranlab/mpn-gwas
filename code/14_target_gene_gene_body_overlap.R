library(tidyverse)
library(data.table)
library(GenomicRanges)
library(BuenColors)
library(annotables)
library(qvalue)

"%ni%" <- Negate("%in%")

# Load protein coding annotations
grch38.pc <- grch38 %>% filter(biotype == "protein_coding")

# Gene body annotations
gencode_gr <- readRDS("../data/annotations/gencode_gene_annotations.rds")

# Read in fine-mapped CS
CS.df<- fread("../data/abf_finemap/MPN_arraycovar_meta_finngen_r4_abf_cojo_95CS.bed")

# Remove regions with high PP coding variant
VEP.PP10 <- fread("../output/VEP/MPN_arraycovar_meta_finngen_r4.coding_variants_PP01.tsv") %>% filter(PP>0.10 | sentinel=="yes")
CS.df <- CS.df %>% filter(region %ni% VEP.PP10$region)

# Merge CS with gene bodies
CS.gr <- GRanges(CS.df); end(CS.gr) <- end(CS.gr)-1
gb_idx <- findOverlaps(gencode_gr,CS.gr)

gb_overlap <- data.frame(
  CS.df[gb_idx@to,c("var","rsid","PP","region","sentinel")],
  gencode_gr[gb_idx@from]
) %>% dplyr::select(-width,-strand) %>% arrange(region,desc(PP))

write.table(gb_overlap,"../output/target_genes/r4_gene_body_target_genes.tsv",
            quote = FALSE, sep = "\t", col.names = T, row.names = F)

# See gene targets of fine-mapped variants
gb_overlap %>% distinct(var,.keep_all = T) %>% .$gene_name %>% unique()
gb_overlap %>% filter(PP>0.1 | sentinel == "yes") %>% distinct(var,.keep_all = T) %>% .$gene_name %>% unique()
