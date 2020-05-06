library(BuenColors)
library(data.table)
library(GenomicRanges)
library(reshape2)
library(ComplexHeatmap)
library(matrixStats)
library(SummarizedExperiment)
library(Matrix)
library(preprocessCore)
library(tidyverse)
library(qvalue)
library(scales)
set.seed(1026)

toplot <- TRUE
# Load approximate BF fine-mapped data
CS.df <- fread("../data/abf_finemap/MPN_arraycovar_meta_finngen_r4_abf_cojo_PP0.001.bed")%>%
  setNames(.,c("seqnames","start","end","PP","region","var","rsid"))

CS.gr <- GRanges(CS.df); end(CS.gr) <- end(CS.gr)-1

# ATAC peaks/counts
peaksdf <- fread("../data/atac/26August2017_EJCsamples_allReads_250bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts.df <-  data.matrix(fread("../data/atac/26August2017_EJCsamples_allReads_250bp.counts.txt"))

# Log2 cpm normalize
cpm <- round(sweep(counts.df, 2, colSums(counts.df), FUN="/") * 1000000, 1)
log2cpm <- log2(cpm+1)
# Min / max scale
log2cpm.minmax <- log2cpm / rowMax(log2cpm)

# Select for myeloid populations
myeloid.counts <- as.data.frame(log2cpm.minmax) %>% 
  dplyr::select("HSC", "MPP", "CMP", "MEP", "Ery","Mega") %>% as.matrix()

# Find overlaps between peaks and CS
pp_threshold = 0.01
CS.df.PP <- CS.df %>% filter(PP > pp_threshold) %>% distinct(var,.keep_all = T)
idx <- findOverlaps(peaks,GRanges(CS.df.PP))
fm.counts <- log2cpm.minmax[idx@from,]

# Percent of variants falling in a peak
nrow(fm.counts) / nrow(CS.df.PP)

# How many variants have max ATAC in a myeloid cell type?
fm.myeloid.counts <- myeloid.counts[idx@from,]
table(rowMax(fm.myeloid.counts) == 1)

# Compare with low PP variants 
if (FALSE){
  pp_threshold = 0.01
  CS.df.lowPP <- CS.df %>% filter(PP < pp_threshold) %>% distinct(var,.keep_all = T)
  idx <- findOverlaps(peaks,GRanges(CS.df.lowPP))
  lowPP.counts <- log2cpm.minmax[idx@from,]
  nrow(lowPP.counts) / nrow(CS.df.lowPP)
  
  chisq.test(matrix(c(nrow(CS.df.PP),nrow(CS.df.lowPP),nrow(fm.counts),nrow(lowPP.counts)),nrow=2))
  fisher.test(matrix(c(nrow(CS.df.PP),nrow(CS.df.lowPP),nrow(fm.counts),nrow(lowPP.counts)), nrow = 2))
}

# Construct mega data.frame with peaks, cpm, and finemap info
mega_df <- data.frame(
  CS.df.PP[idx@to,c("var","rsid","region","PP","pvalue","region_rank")],
  peaks[idx@from],
  log2cpm.minmax[idx@from,]
) %>% dplyr::select(-strand,-width) %>% unique()

# Write all
if (FALSE){
  write.table(mega_df,file="../output/atac_overlap/r4_atac_overlap_all_PP01.tsv",quote = FALSE, sep = "\t", col.names = T, row.names = F)
}

# Annotate with PP bins
range(mega_df$PP)
counts.PP <- as.data.frame(cut(mega_df$PP, c(0.01, 0.05, 0.1, 0.25,0.75, 1.0)))

names(counts.PP) <- "PPbin"

# Cluster by ATAC and plot heatmap
km <- kmeans(fm.counts, centers = 5, nstart = 10000) 
km.cluster <- factor(km$cluster, levels = c("3","2","5","4","1")) 

if (toplot){
  pdf(file=paste0("../output/atac_overlap/r4_MPN_CML.abf.pp",pp_threshold*100,".250bp_counts_allpops.atac.heatmap.pdf"), width = 4, height = 4)
  par(cex.main=0.8,mar=c(1,1,1,1))
  hm <- Heatmap(fm.counts, col=as.character(jdb_palette("brewer_spectra",type="continuous")),
                cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = TRUE,
                row_names_gp = gpar(fontsize = 0),
                column_names_gp = gpar(fontsize = 6),
                split = km.cluster, show_heatmap_legend = T,
                name = "Accessibility")
  ha1 <- rowAnnotation(df = as.data.frame(counts.PP),
                       col = list(PPbin = c("(0.01,0.05]" = jdb_palette("brewer_red")[2],
                                            "(0.05,0.1]" = jdb_palette("brewer_red")[3],
                                            "(0.1,0.25]" = jdb_palette("brewer_red")[5],
                                            "(0.25,0.75]" = jdb_palette("brewer_red")[7],
                                            "(0.75,1]" = jdb_palette("brewer_red")[9])),
                       width = unit(0.8, "cm"), show_legend = T)
  hm + ha1
  
  dev.off()
}

