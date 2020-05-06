library(data.table)
library(tidyverse)
library(BuenColors)
library(matrixStats)
library(ComplexHeatmap)
library(Matrix)
library(SummarizedExperiment)
library(scales)
library(annotables)
"%ni%" <- Negate("%in%")

target_genes <- read.table("../output/target_genes/r4_target_genes_scored_noABC_ties.txt")$V1

# Read in Corces RNA data
rna <- data.frame(fread(paste0("../data/rna/16populations_RNAcounts.txt")))
genes <- rna[,1]
counts <- data.matrix(rna[,-1])
cpm <- sweep(counts, 2, Matrix::colSums(counts), FUN="/") * 1000000
colnames(cpm) <- c("HSC", "MPP", "LMPP", "CLP", "GMP-A", "GMP-B", "GMP-C", "CMP", "MEP", "NK", "CD4", "CD8", "B", "pDC", "Mono", "Ery")
rownames(cpm) <- genes

# Filter genes
keepGenes <- which(rowMeans(counts) > 1)
length(keepGenes)
cpm <- cpm[keepGenes, ]

log2cpm <- log2(cpm+1)

# Load protein coding annotations
grch38.pc <- grch38 %>% filter(biotype == "protein_coding")
log2cpm <- log2cpm[rownames(log2cpm) %in% grch38.pc$symbol,] 

# Min / max scale
log2cpm.minmax <- log2cpm / rowMax(log2cpm)

scaleRows <- function(x) {
  rmeans <- rowMeans(x)
  x <- sweep(x, 1, rmeans)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
}
# z-score the rows (x - rowmeans / sd)
# log2cpm.minmax <- scaleRows(log2cpm)

# Subset to target genes and cell types
if (FALSE){
  ct <- c("HSC","MPP","CMP","LMPP","GMP-A","CLP")
  tg_cpm <- as.data.frame(log2cpm.minmax[rownames(log2cpm.minmax) %in% target_genes,]) %>% 
    dplyr::select(ct) %>% as.matrix()
}

tg_cpm <- as.data.frame(log2cpm.minmax[rownames(log2cpm.minmax) %in% target_genes,]) %>% as.matrix()

# Gap statistic
findBest <- cluster::clusGap(tg_cpm, FUN = kmeans, K.max = 10, B = 100, nstart = 100)
qplot(1:10,findBest$Tab[,3])

# Make final scree plot
df <- data.frame(K = 1:10, Gap = findBest$Tab[,3])
p1 <-ggplot(df, aes(x = K, y = Gap)) + geom_point() +
  pretty_plot(fontsize = 8) + L_border() + labs(x = "K-means", y = "Gap Statistic") +
  geom_vline(xintercept = 5, linetype = 2) +
  scale_x_continuous(breaks= pretty_breaks())
cowplot::ggsave2(p1,file="../output/target_genes/r4_scored_noABC_target_genes_expression_heatmap_gap.pdf",height=2,width=2)

# Cluster by ATAC and plot heatmap
km <- kmeans(tg_cpm, centers = 4, nstart = 10000) 
km.cluster <- factor(km$cluster, levels = c("1","4","2","3"))

pdf(file=paste0("../output/target_genes/r4_scored_noABC_min2_target_genes_expression_heatmap.pdf"), width = 3, height = 2.7)
  Heatmap(tg_cpm, col=as.character(jdb_palette("brewer_spectra",type="continuous")),
          cluster_rows = TRUE, cluster_columns = FALSE, row_title = NULL,
          cluster_row_slices = FALSE, 
          row_names_gp = gpar(fontsize = 6),
          column_names_gp = gpar(fontsize = 6),
          border=TRUE,
          split = km.cluster, 
          show_heatmap_legend = FALSE,
          name = "RNA")
dev.off()
