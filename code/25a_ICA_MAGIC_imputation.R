require(Matrix)
require(Seurat)
require(rhdf5)
require(data.table)
require(future.apply)
library(tidyverse)
library(Rmagic)
library(viridis)
library(phateR)
library(BuenColors)

# Analyze MAGIC results ---------------------------------------------------
ICA_path <- "file path for scRNAseq ICA data"
data <- t(fread(ICA_path))

h5file <- "h5 file path"

gene_names  <-h5read(h5file,"row_attrs/GeneNames")
clusters <- h5read(h5file,"col_attrs/res.0.5")
cluster_ids <- h5read(h5file,"col_attrs/res.0.5_l")
clusters_names <- cluster_ids[clusters]

# Read MPN target genes
mpn <- read.table("/broad/sankaranlab/MPN_GWAS/metaGWAS/target_genes/r4_target_genes_scored_withABC_ties.txt")[,1]  %>% as.character()

# HSC signature
hsc <- c("CD34","HLF","CRHBP")

# For a specified geneset, returns the top Louvain clusters ranked by the expression of the geneset
get_top_clusters <- function(geneset){
  # Retrieve indices with MPN or HSC gene
  idx <- which(rownames(data) %in% geneset)
  names <- rownames(data)[idx]
  
  # Subset matrix to the relevant genes
  subset <- data[idx,]
  rownames(subset) <- names
  
  # z-score each marker gene, then average them
  zscores <- apply(subset, 1, scale)
  meta_z <- apply(zscores ,1,mean)
  # meta_z <- rowSums(zscores) / sqrt(ncol(zscores))
  
  # Compare scores across clusters
  df <- data.frame(id = seq(length(meta_z)),
                   cluster = clusters_names,
                   signature = meta_z)
  
  return(df)
}

hsc_clusters <- get_top_clusters(geneset=hsc) %>% dplyr::rename(hsc_signature="signature")
mpn_clusters <- get_top_clusters(geneset=mpn) %>% dplyr::rename(mpn_signature="signature")
merged <- inner_join(hsc_clusters,mpn_clusters %>% dplyr::select(-cluster),by="id")

cor.test(merged$hsc_signature,merged$mpn_signature,method="spearman")

p1 <- ggplot(merged,aes(x=hsc_signature,y=mpn_signature,color=hsc_signature))+ 
  geom_point_rast(raster.dpi=400,size=0.01)+
  scale_color_viridis(option="B") +
  pretty_plot(fontsize = 8)+ L_border()
ggsave2(p1,file="../output/target_genes/scRNAseq/MPN_HSC_scores_imputed_singlecell_correlation.pdf", height = 2,width =3)

cluster_sizes <- merged %>% group_by(cluster) %>% summarise(n = n())

