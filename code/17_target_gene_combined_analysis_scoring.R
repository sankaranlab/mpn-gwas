library(tidyverse)
library(data.table)
library(GenomicRanges)
library(BuenColors)
library(annotables)
library(qvalue)
library(plotly)
library(Vennerable)
library(preprocessCore)
library(ComplexHeatmap)
library(gtools)

"%ni%" <- Negate("%in%")

# Load protein coding annotations
grch38.pc <- grch38 %>%filter(biotype == "protein_coding")

# Read in target genes from different analyses
gb_overlap <- fread("../output/target_genes/r4_gene_body_target_genes.tsv")
pg_overlap <- fread("../output/target_genes/r4_atac_rna_target_genes.tsv")
myeloid_overlap <- fread("../output/target_genes/r4_pchic_myeloid_target_genes.tsv")

magma <- fread("../output/target_genes/magma/magma.genes.out")
magma_genes <- magma %>% filter(P < 2.633e-6, SYMBOL %in% grch38.pc$symbol) %>% .$SYMBOL

pp_threshold <- 0.01

# Read in fine-mapped CS
CS.df<- fread("../data/abf_finemap/MPN_arraycovar_meta_finngen_r4_abf_cojo_95CS.bed")
CS.gr <- GRanges(CS.df); end(CS.gr) <- end(CS.gr)-1
sentinels <- CS.df %>% filter(sentinel == "yes")
gw_sig_regions <- CS.df %>% filter(sentinel == "yes",pvalue < 5e-8) %>% .$region

# Gene body genes - only keep unique genes
gb <- gb_overlap %>% filter(PP> pp_threshold | sentinel == "yes") %>% dplyr::rename("gene"=gene_name) %>%
  dplyr::select(var,rsid,PP,region,sentinel,gene) %>% 
  arrange(desc(PP)) %>% 
  distinct(gene,.keep_all = T)  %>%
  mutate(type = "genebody")

# ATAC-RNA correlation genes
pg <- pg_overlap %>% filter(PP> pp_threshold| sentinel == "yes") %>% 
  dplyr::select(var,rsid,PP,region,sentinel,gene) %>% 
  arrange(desc(PP)) %>%
  distinct(gene,.keep_all = T)  %>%
  mutate(type = "atac_rna")

# PCHiC genes
myeloid <- myeloid_overlap %>% filter(PP> pp_threshold| sentinel == "yes") %>% 
  dplyr::select(var,rsid,PP,region,sentinel,gene) %>% 
  arrange(desc(PP)) %>%
  distinct(gene,.keep_all = T)  %>%
  mutate(type = "myeloid_pchic")

# Add in coding genes
VEP.PP10 <- fread("../output/VEP/MPN_arraycovar_meta_finngen_r4.coding_variants_PP01.tsv") %>% 
  filter(PP > pp_threshold| sentinel =="yes") %>% 
  dplyr::rename(gene="SYMBOL") %>%  dplyr::select(var,rsid,PP,region,sentinel,gene) %>% 
  arrange(desc(PP)) %>%
  distinct(gene,.keep_all = T)  %>%
  mutate(type = "coding") 

# Add coding genes
coding_genes <- VEP.PP10 %>% filter(PP>0.10) 

# Combine and filter for protein-coding genes
combined_variants <- bind_rows(gb,pg,myeloid,coding_genes) %>% unique() %>% 
  filter(gene %in% grch38.pc$symbol)

# Remove genes from HLA region
hla_bounds <- c(28866528,33775446)
hla_exclude <- grch38 %>% filter(chr==6  & start > hla_bounds[1] & start < hla_bounds[2]) %>% .$symbol
combined_variants <- combined_variants %>% filter(gene %ni% hla_exclude)

# Add MAGMA significant genes
magma.df <- data.frame(gene = magma_genes, type = "magma")
magma.df <- left_join(magma.df,combined_variants[,c("region","gene")],by="gene") %>% unique()

# Incorporate Grinfeld somatic gene information
somatic_genes <- read.table("../data/annotations/somatic_mutated_genes_MPN_Grinfeld.txt")[,1]
somatic.df <- left_join(data.frame(gene = somatic_genes,type="somatic"),
                        combined_variants[,c("region","gene")],by="gene") %>% unique()

combined_genes <- bind_rows(combined_variants %>% dplyr::select(region,gene,type) %>% unique(),
                            magma.df,
                            somatic.df) %>% drop_na()

# Score
scores <- combined_genes %>% group_by(gene,region) %>% summarise(score = n()) %>% as.data.frame() %>% 
  arrange(region,desc(score))
# coding variant counts for 10 points
scores <- scores %>% mutate(score = ifelse(gene %in% coding_genes$gene,score + 9, score))

# Pick top scoring genes from each locus. Multiple genes are included for ties
scores_topregion <- scores %>% group_by(region) %>% filter(score == max(score)) %>%
  ungroup() %>% unique() %>% as.data.frame() %>%
  filter(score > 1) 

if (TRUE){
  # Add chromosome cytoband information
  cytobands.gr <- fread("../data/annotations/cytoBand.txt") %>% 
    setNames(.,c("seqnames","start","end","locus","stain")) %>% mutate(locus = paste0(gsub("chr","",seqnames),locus)) %>% 
    GRanges()
  
  scores_locus <- merge(scores_topregion,grch37[,c("symbol","chr","start","end")],by.x="gene",by.y="symbol") %>%
    filter(chr %in% seq(1,22)) %>%
    mutate(chr = paste0("chr",chr),end=start) %>% unique()
  scores_locus.gr <- GRanges(scores_locus)
  locus_idx <- findOverlaps(scores_locus.gr,cytobands.gr)
  scores_locus$locus <- cytobands.gr[locus_idx@to]$locus
  scores_topregion <- merge(scores_topregion,scores_locus[,c("region","locus")]) %>% unique()
  
  scores_topregion<- scores_topregion %>% distinct(region,gene,.keep_all = TRUE)
}

fwrite(data.frame(gene=unique(scores_topregion$gene)),sep="\t",file="r4_target_genes_scored_withABC_ties.txt",col.names=FALSE)

# Make gene x evidence heatmap
# Assign one locus to each gene
scores_topregion_sumstats <- merge(scores_topregion,sentinels[,c("rsid","region","pvalue")],by="region") %>%
  arrange(region) %>% unique()
scores_topregion_sumstats <- scores_topregion_sumstats %>% group_by(gene) %>% filter(score == max(score)) %>% ungroup() %>% unique()

combined_genes_score <- combined_genes[,c("gene","type")] %>% 
  mutate(score = ifelse(type == "coding",2,1)) %>% unique()
gene_matrix <- data.matrix(pivot_wider(combined_genes_score,names_from=type,values_from=score)[,-1])
rownames(gene_matrix) <- pivot_wider(combined_genes_score,names_from=type,values_from=score)$gene
gene_matrix[is.na(gene_matrix)] <- 0

if (FALSE){
  # All genes
  cor(gene_matrix)
  Heatmap(gene_matrix, col=c("white","firebrick","dodgerblue"),
          cluster_rows = FALSE, cluster_columns = FALSE, 
          show_column_names = TRUE,
          row_names_gp = gpar(fontsize = 6),
          column_names_gp = gpar(fontsize = 6),
          show_heatmap_legend = FALSE)
}

# Subset to selected genes
gene_matrix_filtered <- gene_matrix[unique(scores_topregion_sumstats$gene),]

totscore <- rowSums(gene_matrix_filtered)
totscore[coding_genes$gene] <-  totscore[coding_genes$gene] - 1 

ha = rowAnnotation(
  Locus = anno_text(scores_topregion_sumstats$locus, location = 0.5,just="center",
                    gp = gpar(border="black",fontsize=6)),
  width = max_text_width(scores_topregion_sumstats$locus)*1.2,
  total_score = anno_simple(totscore,
                            pch = c(as.character(totscore)),
                            col = c("1" = jdb_palette("brewer_red")[1],
                                    "2" = jdb_palette("brewer_red")[3],
                                    "3" = jdb_palette("brewer_red")[5],
                                    "4" = jdb_palette("brewer_red")[7]),
                            gp = gpar(col = "black"),
                            pt_gp = gpar(fontsize=6)),
  annotation_name_gp = gpar(fontsize=6), annotation_name_side = "top")


pdf(file=paste0("../output/target_genes/r4_target_gene_scored_noABC_min2_annotation_heatmap.pdf"), width = 3, height = 4)
par(cex.main=0.8,mar=c(1,1,1,1))
  Heatmap(gene_matrix_filtered, col=c("white","firebrick3","deepskyblue3"),
          column_order = c("genebody","magma","atac_rna","myeloid_pchic","somatic","coding"),
          cluster_rows = FALSE, cluster_columns = FALSE, 
          show_column_names = TRUE,
          row_names_gp = gpar(fontsize = 6),
          column_names_gp = gpar(fontsize = 6),
          column_names_side = "top",
          row_names_side = "left",
          border=TRUE,rect_gp = gpar(col = "black", lwd = 1),
          show_heatmap_legend = FALSE,
          left_annotation = ha)
dev.off()

table(combined_variants$region) %>% length()
gw_sig_combined_variants <- combined_variants %>% filter(region %in% gw_sig_regions)
