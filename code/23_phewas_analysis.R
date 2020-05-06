library(data.table)
library(tidyverse)
library(qvalue)
library(Matrix)
library(ComplexHeatmap)
library(SummarizedExperiment)
library(BuenColors)
library(annotables)
library(ggrastr)
library(ggrepel)
"%ni%" <- Negate("%in%")

# Write list of variants for pheWAS input ---------------------------------
if (FALSE){
  # Add leading 0 to chromosomes under 10
  CS.df <- fread("../data/abf_finemap/MPN_arraycovar_meta_finngen_r4_abf_cojo_PP0.001_annotated.bed") %>% 
    mutate(chr = as.integer(gsub("chr","",seqnames))) %>% mutate(var = ifelse(chr < 10,paste0("0",var),var)) %>%
    dplyr::select(var) %>% unique()
  fwrite(CS.df,file="../output/phewas/r4_MPN_PP001_phewas_input_leading_zeros.txt",sep="\t",col.names = F)
  
  # Non-leading 0
  CS.df <- fread("../data/abf_finemap/MPN_arraycovar_meta_finngen_r4_abf_cojo_PP0.001_annotated.bed")  %>%
    dplyr::select(var) %>% unique()
  fwrite(CS.df,file="../output/phewas/r4_MPN_PP001_phewas_input_nonleading_zeros.txt",sep="\t",col.names = F)
  
}

# SAIGE pheWAS ------------------------------------------------------------
# Read phewas results
icd <- fread("zcat < ../output/phewas/SAIGE.allphenos.mpn_r4.txt.gz") %>%
  mutate(var = paste(CHROM,POS,REF,ALT,sep=":"), maf = ifelse(af < 0.5, af, 1-af),
         expected_case_minor_AC =  2 * maf * num_cases, logp = -log10(pval),
         var = paste0(CHROM,":",POS,"_",REF,"_",ALT)) %>%
  filter(expected_case_minor_AC > 50)

# Read and merge pheno codes 
phecodes <- readxl::read_xlsx("../output/phewas/saige-phenotype-information.xlsx") %>% 
  dplyr::rename(pheno = `Phenotype Description`,category=`Phenotype Category`) %>%
  filter(category %ni% c("pregnancy complications","injuries & poisonings","sense organs","symptoms"))
phecodes$PheCode <- as.double(phecodes$PheCode)
icd <- merge(icd,phecodes[,c("PheCode","pheno","category")],by.x="phecode",by.y="PheCode")

# Combine with sum stats 
CS.df <- fread("../data/abf_finemap/MPN_arraycovar_meta_finngen_r4_abf_cojo_PP0.001_annotated.bed")
merged <- merge(icd,CS.df[,c("var","rsid","PP","pvalue","region","region_rank","effect","stderr")],by="var") %>%
  mutate(CHROM=paste0("chr",CHROM))

# Remove MPNs and polycythemia vera due to redundancy
merged <- merged %>% filter(pheno %ni% c("Myeloproliferative disease","Polycythemia vera","Polycythemia, secondary"))

# Exclude HLA region
hla_bounds <- c(28866528,33775446)
toexclude <- merged %>% filter(CHROM=="chr6", POS > hla_bounds[1],POS < hla_bounds[2]) %>% .$var
merged <- merged %>% filter(var %ni% toexclude)

# Add nearest gene
# Read in gene body annotations
gencode_gr <- readRDS("../../bcx-finemap/data/annotations/gencode_filtered.rds")

# Find nearest gene
merged.gr <- makeGRangesFromDataFrame(merged,keep.extra.columns = T,seqnames.field="CHROM",
                                      start.field = "POS",end.field = "POS")
merged$nearest_gene <- gencode_gr[nearest(merged.gr,gencode_gr,ignore.strand=TRUE),]$gene_name

# Set bonferroni threshold
highPP <- merged %>% filter(PP > 0.10| region_rank ==1 )
num_variants <- highPP %>% .$var %>% unique %>% length() 
p_threshold <- 0.05 / (length(unique(merged$pheno)))

# Heatmap
if (FALSE){
  # Read in cytoband information
  cytobands.gr <- fread("../data/annotations/cytoBand.txt") %>% 
    setNames(.,c("seqnames","start","end","locus","stain")) %>% mutate(locus = paste0(gsub("chr","",seqnames),locus)) %>% 
    GRanges()
  
  # Orient by risk allele
  sumstats <- "file path to summary statistics"
  risk_oriented <- fread(sumstats)
  risk_icd <- inner_join(icd,risk_oriented[,c("MarkerName","risk","nonrisk")],by=c("var"="MarkerName"))
  flipped <- risk_icd %>% filter(ALT != risk) %>% mutate(beta = -1*beta)
  risk_merged <- risk_icd %>% filter(ALT == risk) %>% bind_rows(.,flipped) %>% mutate(Z = beta/sebeta)
  
  # Group by variants, take top variant per region
  selected <-   sig_vars <- merged %>% filter(pval < p_threshold) %>% 
    filter(PP > 0.10| region_rank ==1 ) %>% 
    group_by(region,pheno) %>% dplyr::slice(which.max(PP)) %>% ungroup() %>% as.data.frame() %>% arrange(region)
  sig_phenos <-selected %>% .$pheno %>% unique()
  sig_vars <- selected %>% .$rsid %>% unique()
  gene_annot <- selected %>% distinct(rsid,.keep_all = TRUE) %>% .$nearest_gene

  # Bin the p-values
  # risk_merged$bin <- cut(risk_merged$pval, c(1, 0.01, p_threshold, 5e-8, 1e-20,0))
  # Bin z-scores
  risk_merged$bin <- cut(risk_merged$Z, c(min(risk_merged$Z)-0.01,-2.58 ,2.58, 4.084,8,max(risk_merged$Z)+0.01))
  table(risk_merged$bin)
  2*pnorm(-abs(4.084))
  
  logp_matrix <- risk_merged %>% filter(pheno %in% sig_phenos) %>% dplyr::select(ID,pheno,bin) %>%
    pivot_wider(names_from=pheno,values_from=bin) %>% as.data.frame()
  rownames(logp_matrix) <- logp_matrix[,1]
  logp_matrix <- logp_matrix[,-1] %>% as.matrix()
  logp_matrix <- logp_matrix[sig_vars,]
  logp_matrix[is.na(logp_matrix)] <-levels(risk_merged$bin)[length(levels(risk_merged$bin))]
  
  # Get cytobands
  CS.subset <- CS.df %>% filter(rsid %in% rownames(logp_matrix)) %>% arrange(region) %>% GRanges()
  locus_idx <- findOverlaps(CS.subset,cytobands.gr)
  CS.subset$locus <- cytobands.gr[locus_idx@to]$locus

  logp_matrix <- logp_matrix[CS.subset$rsid,]
  
  # Get pheno categories
  categories <- icd %>% distinct(pheno, .keep_all = T) %>%
    filter(pheno %in% colnames(logp_matrix)) %>% .$category
  table(categories)

  # Heatmap
  ha = rowAnnotation(
    Locus = anno_text(CS.subset$locus, location = 0.5,just="center",
                      gp = gpar(border="black",fontsize=6)),
    Nearest_gene = anno_text(gene_annot, location = 0.5,just="center",
                             gp = gpar(border="black",fontsize=6)),
    width = max_text_width(CS.subset$locus),
    annotation_name_gp = gpar(fontsize=6), annotation_name_side = "top")
  
  block_anno = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = jdb_palette("lawhoops")),
                                                  labels = unique(categories), 
                                                  labels_rot = 90,
                                                  labels_gp = gpar(col = "white", fontsize = 8)),
                                 height = max_text_width(categories))
  
  pdf(file=paste0("../output/phewas/r4_phewas_heatmap_zscores.pdf"), width = 8, height = 4.5)
  
  colors = structure(jdb_palette("solar_flare")[c(9,7,6,5,3)], 
                     names =rev(levels(risk_merged$bin)))
  # colors = jdb_palette("solar_extra")
  Heatmap(logp_matrix, col=colors,
          cluster_rows = FALSE, cluster_columns = FALSE, 
          show_column_names = TRUE, column_title = NULL,
          row_names_gp = gpar(fontsize = 6),
          column_names_gp = gpar(fontsize = 6),
          column_names_rot = 45,
          top_annotation = block_anno,
          column_split = factor(categories,levels = unique(categories)),
          show_heatmap_legend = T,
          border=TRUE,rect_gp = gpar(col = "black", lwd = 1),
          name = "z-score",
          row_names_side = "left",
          left_annotation = ha)
  dev.off()
}  
