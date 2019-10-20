library(tidyverse)
library(data.table)
library(GenomicRanges)
library(BuenColors)
library(annotables)
library(qvalue)
library(plotly)
library(Vennerable)
library(preprocessCore)

"%ni%" <- Negate("%in%")

# Load protein coding annotations
grch38.pc <- grch38 %>%filter(biotype == "protein_coding")

# Read in target genes from different analyses
gb_overlap <- fread("../output/target_genes/gene_body_target_genes.tsv")
pg_overlap <- fread("../output/target_genes/atac_rna_target_genes.tsv")
myeloid_overlap <- fread("../output/target_genes/pchic_myeloid_target_genes.tsv")
cd34_overlap <- fread("../output/target_genes/pchic_cd34_target_genes.tsv")

pp_threshold <- 0.1

# Read in fine-mapped CS
CS.df<- fread("../data/abf_finemap/MPN_CML_abf_cojo_95CS.bed")
CS.gr <- GRanges(CS.df); end(CS.gr) <- end(CS.gr)-1
CS.df %>% filter(sentinel == "yes",pvalue < 5e-8) %>% .$region -> gw_sig_regions

# Gene body genes
gb <- gb_overlap %>% filter(PP> pp_threshold | sentinel == "yes") %>% dplyr::rename("gene"=gene_name) %>% dplyr::select(var,rsid,PP,region,sentinel,gene) %>% distinct(var,gene,.keep_all = T) 

# ATAC-RNA correlation genes
pg <- pg_overlap %>% filter(PP> pp_threshold | sentinel == "yes") %>% group_by(var) %>%
  filter(pvalue == min(pvalue)) %>% dplyr::select(var,rsid,PP,region,sentinel,gene)

# PCHiC genes
myeloid <- myeloid_overlap %>% filter(PP> pp_threshold | sentinel == "yes") %>% group_by(var) %>%
  filter(maxscore == max(maxscore)) %>% dplyr::select(var,rsid,PP,region,sentinel,gene) %>% distinct(var,gene,.keep_all = T)
cd34 <- cd34_overlap %>% filter(PP> pp_threshold | sentinel == "yes")  %>% group_by(var) %>%
  filter(CD34 == max(CD34))   %>% dplyr::select(var,rsid,PP,region,sentinel,gene)%>% distinct(var,gene,.keep_all = T)

# Add in coding genes
VEP.PP10 <- fread("../output/VEP/coding_variants_PP01.tsv") %>% filter(PP>pp_threshold| sentinel =="yes") %>% dplyr::rename(gene="SYMBOL") %>%  dplyr::select(var,rsid,PP,region,sentinel,gene) 

# Combine and filter for protein-coding genes
combined_variants <- bind_rows(gb,pg,myeloid,cd34,VEP.PP10) %>% unique() %>% 
  filter(gene %in% grch38.pc$symbol)
combined_genes <- unique(combined_variants$gene)

# Remove genes from HLA region
hla_bounds <- c(28866528,33775446)
hla_exclude <- grch38 %>% filter(chr==6  & start > hla_bounds[1] & start < hla_bounds[2]) %>% .$symbol

combined_genes <- combined_genes[combined_genes %ni% hla_exclude]
combined_variants <- combined_variants %>% filter(gene %ni% hla_exclude)
gw_sig_combined_variants <- combined_variants %>% filter(region %in% gw_sig_regions)

# Write list of target genes
write.table(combined_genes,file="../output/target_genes/target_genes_combined.PP10.txt",quote = FALSE, sep = "\t", col.names = F, row.names = F)

# Plot Venn diagram of overlap
gene_venn <- Venn(list(gene_body=gb$gene,peak_gene=pg$gene,PCHiC=unique(c(myeloid$gene,cd34$gene))))  
pdf("../output/target_genes/targetgene_venn.pdf")
plot(gene_venn,doWeights=TRUE)
dev.off()


# Create master table of variant to gene ----------------------------------

# Merge variants with sumstats and ATAC
CS.df <- fread("../data/abf_finemap/MPN_CML_abf_cojo_95CS.bed")
merged <- merge(combined_variants,CS.df[,c("var","pvalue","maf","effect","stderr")],by="var") %>% arrange(region,desc(PP))

# Peaks data
peaksdf <- fread("../data/atac/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
# Counts data
counts.df <-  data.matrix(fread("../data/atac/29August2017_EJCsamples_allReads_500bp.counts.txt"))

# Subset to only good peaks and counts
if (F){
  n = 0.75
  # Quantile normalize the counts matrix
  counts.df[1:dim(counts.df)[1],1:dim(counts.df)[2]] <- round(normalize.quantiles(as.matrix(counts.df)),2)
  keep <- apply(counts.df,1,max) > mean(apply(counts.df,2,function(x) {quantile(x,n)}))
  counts.df <- counts.df[keep,]
  
  peaks <- peaks[keep,]
}
# Log2 cpm normalize
cpm <- round(sweep(counts.df, 2, colSums(counts.df), FUN="/") * 1000000, 1)
log2cpm <- log2(cpm+1)

# Select for myeloid populations
myeloid.counts <- as.data.frame(log2cpm) %>% 
  dplyr::select("HSC", "MPP", "CMP", "MEP", "Ery","Mega") %>% as.matrix()


merged$seqnames <-paste0("chr",gsub("_.*","",str_split_fixed(merged$var,":",2))[,1]) 
merged$end <- merged$start <-gsub("_.*","",str_split_fixed(merged$var,":",2))[,2]
idx <- findOverlaps(peaks,GRanges(merged))

# Construct mega data.frame with peaks, cpm, and finemap info
mega_df <- data.frame(
  var = merged[idx@to,"var"],
  myeloid.counts[idx@from,]
)

all_merged <- left_join(merged,mega_df,by="var") %>% dplyr::select(-start,-end) %>% unique()
write.table(all_merged,file="../output/target_genes/target_genes_variant_list.tsv",
            quote = FALSE, sep = "\t", col.names = T, row.names = F)
