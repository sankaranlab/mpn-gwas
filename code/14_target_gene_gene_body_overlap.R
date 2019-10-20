library(tidyverse)
library(data.table)
library(GenomicRanges)
library(BuenColors)
library(annotables)
library(qvalue)

"%ni%" <- Negate("%in%")

# Load protein coding annotations
# grch38.pc <- grch38 %>% filter(biotype == "protein_coding" | biotype == "antisense_RNA")
grch38.pc <- grch38 %>% filter(biotype == "protein_coding")

if (FALSE){
  # Read in gene body annotations
  gencode_gtf <- fread(paste0("zcat <","../data/annotations/gencode.v28lift37.basic.annotation.gtf.gz"))
  colnames(gencode_gtf) <- c("chr","annotation_source","feature_type","start","end","score","strand","phase","additional_info")
  
  # Only keep gene name, remove special characters
  gencode_gtf$gene_name <- str_replace_all(gsub(";.*","",gsub(".*gene_name ","",gencode_gtf$additional_info)), "[[:punct:]]", "")
  gencode_gtf$gene_id <- str_replace_all(gsub("\\..*","",gsub(".*gene_id ","",gencode_gtf$additional_info)), "[[:punct:]]", "")
  
  # Remove RP genes, filter for only protein coding genes
  gencode_gtf_filtered <- gencode_gtf %>% filter(!str_detect(gene_name, "^RP")) %>% filter(gene_id %in% grch38.pc$ensgene) 
  
  # Convert to GRanges
  gencode_gr <- makeGRangesFromDataFrame(gencode_gtf_filtered,seqnames = "chr", start.field = "start", end.field = "end",keep.extra.columns = TRUE)
  
  saveRDS(gencode_gr,file="../data/annotations/gencode_gene_annotations.rds")
}
gencode_gr <- readRDS("../data/annotations/gencode_gene_annotations.rds")

# Read in fine-mapped CS
CS.df<- fread("../data/abf_finemap/MPN_CML_abf_cojo_95CS.bed")

# Remove regions with high PP coding variant
VEP.PP10 <- fread("../output/VEP/coding_variants_PP01.tsv") %>% filter(PP>0.10 | sentinel=="yes")
CS.df <- CS.df %>% filter(region %ni% VEP.PP10$region)

# Merge CS with gene bodies
CS.gr <- GRanges(CS.df); end(CS.gr) <- end(CS.gr)-1
gb_idx <- findOverlaps(gencode_gr,CS.gr)

gb_overlap <- data.frame(
  CS.df[gb_idx@to,c("var","rsid","PP","region","sentinel")],
  gencode_gr[gb_idx@from]
) %>% dplyr::select(-width,-strand) %>% arrange(region,desc(PP))

write.table(gb_overlap,"../output/target_genes/gene_body_target_genes.tsv",
            quote = FALSE, sep = "\t", col.names = T, row.names = F)

# See gene targets of fine-mapped variants
gb_overlap %>% filter(PP>0.1 | sentinel == "yes") %>% distinct(var,.keep_all = T) %>% .$gene_name %>% unique()
gb_overlap %>% filter(PP>0.25 | sentinel == "yes") %>% distinct(var,.keep_all = T) %>% .$gene_name %>% unique()

