library(tidyverse)
library(data.table)
library(GenomicRanges)
library(BuenColors)
library(annotables)

"%ni%" <- Negate("%in%")

# Load protein coding annotations
grch38.pc <- grch38 %>%filter(biotype == "protein_coding")

# Blood PCHIC
pchic_blood <- fread("zcat < /Volumes/broad_sankaranlab/pCHIC/pCHiC_Blood.tsv.gz")
all_cells <- colnames(pchic_blood)[7:23]
myeloid_cells <- c("Mon","Neu","MK","Ery","Mac0","Mac1","Mac2","EP")

# Only interactions with score ≥ 5 in at least 1 cell type
allcells <- TRUE
if (allcells){
  pchic_myeloid <- pchic_blood %>% dplyr::rename(gene="baitName") %>%
    dplyr::select(oeChr,oeStart,oeEnd,all_cells,gene)
} else{
  pchic_myeloid <- pchic_blood %>% dplyr::rename(gene="baitName") %>%
    dplyr::select(oeChr,oeStart,oeEnd,myeloid_cells,gene)
}

pchic_myeloid <- pchic_myeloid[apply(pchic_myeloid[,4:(ncol(pchic_myeloid)-1)],1,FUN=max) > 5,]
colnames(pchic_myeloid) <- gsub("oe","",colnames(pchic_myeloid))
pchic_myeloid <- pchic_myeloid[pchic_myeloid$gene %in% grch38.pc$symbol,]
pchic_myeloid$Chr <- paste0("chr",pchic_myeloid$Chr)
pchic_myeloid.gr <- GRanges(pchic_myeloid)

# CD34 PCHIC
pchic_cd34 <- readRDS("../data/pchic/cd34_loops.rds")
pchic_cd34.gr <- GRanges(pchic_cd34)
pchic_cd34.gr <- pchic_cd34.gr[pchic_cd34.gr$gene %in% grch38.pc$symbol,]


# Read in fine-mapped CS
CS.df<- fread("../data/abf_finemap/MPN_CML_abf_cojo_95CS.bed")
# Remove regions with high PP coding variant
VEP.PP10 <- fread("../output/VEP/coding_variants_PP01.tsv") %>% filter(PP>0.01 | sentinel == "yes")
CS.df <- CS.df %>% filter(region %ni% VEP.PP10$region)

CS.gr <- GRanges(CS.df); end(CS.gr) <- end(CS.gr)-1

# Merge CS with pchic
myeloid_idx <- findOverlaps(pchic_myeloid.gr,CS.gr)
cd34_idx <- findOverlaps(pchic_cd34.gr,CS.gr)

myeloid_overlap <- data.frame(
  CS.df[myeloid_idx@to,c("var","rsid","PP","region","sentinel")],
  pchic_myeloid.gr[myeloid_idx@from]
) %>% dplyr::select(-width,-strand) %>% arrange(region,desc(PP))
myeloid_overlap$maxscore <- apply(myeloid_overlap[,myeloid_cells],1,FUN=max)

cd34_overlap <- data.frame(
  CS.df[cd34_idx@to,c("var","rsid","PP","region","sentinel")],
  pchic_cd34.gr[cd34_idx@from]
) %>% dplyr::rename(CD34="score") %>% dplyr::select(-width,-strand)  %>% arrange(region,desc(PP))

# See gene targets of fine-mapped variants
myeloid_overlap %>% filter(PP>0.1 | sentinel == "yes")  %>% .$gene %>% unique()

cd34_overlap %>% filter(PP>0.1 | sentinel == "yes") %>% group_by(var) %>%
  filter(CD34 == max(CD34)) %>% .$gene %>% unique()

# Write target gene tables
write.table(myeloid_overlap,"../output/target_genes/pchic_myeloid_target_genes.tsv",
            quote = FALSE, sep = "\t", col.names = T, row.names = F)
write.table(cd34_overlap,"../output/target_genes/pchic_cd34_target_genes.tsv",
            quote = FALSE, sep = "\t", col.names = T, row.names = F)

