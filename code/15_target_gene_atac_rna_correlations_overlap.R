library(tidyverse)
library(data.table)
library(GenomicRanges)
library(BuenColors)
library(annotables)
library(qvalue)

# Load protein coding annotations
grch38.pc <- grch38 %>%
  dplyr::filter(biotype == "protein_coding")

# Read in enhancer gene correlations
pg.df <- fread(paste0("zcat < ", "../data/rna/peakGeneCorrelation.tsv.gz"))
names(pg.df) <- c("chrom","start","end","gene","cor","pvalue")
pg.df$qvalue <- qvalue(pg.df$pvalue)$qvalues
pg.df <- pg.df %>% filter(qvalue < 0.001) %>% filter(gene %in% grch38.pc$symbol)
pg.gr <- GRanges(pg.df)

# Read in fine-mapped CS
CS.df<- fread("../data/abf_finemap/MPN_CML_abf_cojo_95CS.bed")
# Remove regions with high PP coding variant
VEP.PP10 <- fread("../output/VEP/coding_variants_PP01.tsv") %>% filter(PP>0.10 | sentinel == "yes")
CS.df <- CS.df %>% filter(region %ni% VEP.PP10$region)

CS.gr <- GRanges(CS.df); end(CS.gr) <- end(CS.gr)-1

# Merge CS with ATAC-RNA correlations
pg_idx <- findOverlaps(pg.gr,CS.gr)

pg_overlap <- data.frame(
  CS.df[pg_idx@to,c("var","rsid","PP","region","sentinel")],
  pg.gr[pg_idx@from]
) %>% dplyr::select(-width,-strand) %>% arrange(region,desc(PP))

# See gene targets of fine-mapped variants
pg_overlap %>% filter(PP>0.01 | sentinel == "yes") %>% group_by(var) %>%
  filter(pvalue == min(pvalue)) %>% .$gene %>% unique()

# Write atac_rna gene targets table
write.table(pg_overlap,"../output/target_genes/atac_rna_target_genes.tsv",
            quote = FALSE, sep = "\t", col.names = T, row.names = F)

