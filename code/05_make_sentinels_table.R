library(data.table)
library(tidyverse)
library(GenomicRanges)
library(regioneR)

# Generate gene body annotations
if (FALSE){
  # Read in gene body annotations
  gencode_gtf <- fread(paste0("zcat <","../data/annotations/gencode.v33lift37.basic.annotation.gtf.gz"))
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

# Read in gene body annotations
gencode_gr <- readRDS("../../bcx-finemap/data/annotations/gencode_filtered.rds")

# Read in sentinels
CS.df <-  fread("../data/abf_finemap/MPN_arraycovar_meta_finngen_r4_abf_cojo_PP0.001_annotated.bed") 
sentinels <- CS.df %>% filter(sentinel=="yes") %>%  dplyr::rename(SNP="rsid")

# Merge with summary statistic p-values
sumstats <- "file path to MPN summary statistics"
metal <- fread(sumstats)
sumstats_merged <- merge(sentinels, metal[,c("RSID","RAF","risk","nonrisk","Effect","StdErr")],by.x="SNP",by.y="RSID")
# OR and 95% CI
sumstats_merged$OR <- round(exp(sumstats_merged$Effect),2)
sumstats_merged$OR_95CI <- paste0(round(exp(sumstats_merged$Effect - 1.96*sumstats_merged$StdErr),2),"-",
                                  round(exp(sumstats_merged$Effect + 1.96*sumstats_merged$StdErr),2))

# Nearest gene
# Find nearest gene
sumstats.gr <- GRanges(sumstats_merged)
sumstats_merged$nearest_gene <- gencode_gr[nearest(sumstats.gr,gencode_gr,ignore.strand=TRUE),]$gene_name

# Find genes within 25 kb
window <- 25000
sumstats.expanded.gr <- extendRegions(sumstats.gr,extend.start = window,extend.end=window)
gene_idx <- findOverlaps(sumstats.expanded.gr,gencode_gr)
other_genes <- data.frame(sumstats_merged[gene_idx@from,],
           other_genes = gencode_gr[gene_idx@to]$gene_name) %>% distinct(SNP,other_genes,.keep_all = TRUE) %>%
  mutate(other_genes = as.character(other_genes))
other_genes <- other_genes %>% filter(nearest_gene != other_genes) %>% group_by(SNP) %>% summarize(other_genes = paste(other_genes, collapse = ','))
sumstats_allgenes_merged <- left_join(sumstats_merged,other_genes[,c("SNP","other_genes")],by="SNP") %>% 
  mutate(other_genes = ifelse(nearest_gene == other_genes,"",other_genes)) 
sumstats_allgenes_merged[is.na(sumstats_allgenes_merged)] <- ""

# Add chromosome cytoband information
cytobands.gr <- fread("../data/annotations/cytoBand.txt") %>% 
  setNames(.,c("seqnames","start","end","locus","stain")) %>% mutate(locus = paste0(gsub("chr","",seqnames),locus)) %>% 
  GRanges()

locus_idx <- findOverlaps(sumstats.gr,cytobands.gr)
sumstats_allgenes_merged$locus <- cytobands.gr[locus_idx@to]$locus

# Merge with replication
replication <- fread("../output/replication_joint_statistics/MVP.mpn_expanded.hasJAK2.tsv")
sumstats_allgenes_merged <- left_join(sumstats_allgenes_merged,replication[,c("rsid","meta_pval")],by=c("SNP"="rsid"))

# Output columns
output <- sumstats_allgenes_merged %>% dplyr::rename(Chr = "seqnames",bp="start") %>%arrange(region) %>%
  dplyr::select(locus,SNP,Chr,bp,risk,nonrisk,RAF,OR,OR_95CI,pvalue,meta_pval,nearest_gene,other_genes) 

write.table(output,file="../output/significant_loci_tables/MPN_r4_suggestive_loci.tsv",col=T,row=F,quo=F,sep='\t')
write.table(output %>% filter(meta_pval < 5e-8),file="../output/significant_loci_tables/MPN_r4_genomewide_loci.tsv",col=T,row=F,quo=F,sep='\t')
