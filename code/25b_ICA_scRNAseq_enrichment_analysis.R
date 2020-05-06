require(Matrix)
require(Seurat)
require(rhdf5)
require(data.table)
require(future.apply)
library(tidyverse)

# Run on Unix cluster with 50+ gb memory --------------------------------
if (FALSE){
  setwd("ICA directory")

  print("Reading ICA (bone marrow) dataset")
  raw_h5file <- "ica_bone_marrow_h5.h5"
  h5file <- "ICA_BM.h5"
  
  h5ls(h5file)
  
  data        <-h5read(h5file,"matrix")
  gene_names  <-h5read(h5file,"row_attrs/GeneNames")
  clusters <- h5read(h5file,"col_attrs/res.0.5")
  cluster_ids <- h5read(h5file,"col_attrs/res.0.5_l")
  
  # Read MPN target genes
  mpn <- read.table("../output/target_genes/r4_target_genes_scored_noABC_ties.txt")[,1] %>% as.character()
  
  # HSC signature
  hsc <- c("CD34","HLF","CRHBP")
  
  # Filter for cells that express at least 1 marker gene
  master_subset <- data[which(gene_names %in% c(mpn,hsc)),]
  cells_idx <- which(colSums(master_subset)>0)
  
  # For a specified geneset, returns the top Louvain clusters ranked by the expression of the geneset
  get_top_clusters <- function(geneset){
    # Retrieve indices with MPN or HSC gene
    idx <- which(gene_names %in% geneset)
    names <- gene_names[idx]
    
    # Subset giant matrix to just the relevant genes
    subset <- data[idx,]
    rownames(subset) <- names
    
    # z-score each marker gene, then average them
    zscores <- apply(subset, 1, scale)
    meta_z <- apply(zscores ,1,mean)
    # meta_z <- rowSums(zscores) / sqrt(ncol(zscores))
    
    # Compare scores across clusters
    clusters_names <- cluster_ids[clusters]
    df <- data.frame(id = seq(length(meta_z)),
                     cluster = clusters_names,
                     signature = meta_z)
    
    return(df)
    if (FALSE){
      ranks <- df %>% group_by(cluster) %>% summarise(signature = mean(signature)) %>% arrange(desc(signature))
      return(ranks)
    }
  }
  
  hsc_clusters <- get_top_clusters(geneset=hsc) %>% dplyr::rename(hsc_signature="signature")
  mpn_clusters <- get_top_clusters(geneset=mpn) %>% dplyr::rename(mpn_signature="signature")
  
  # Rank correlation between different gene sets
  signature_file <- "mpn_hsc_signatures file"
  merged <- inner_join(hsc_clusters,mpn_clusters %>% dplyr::select(-cluster),by="id")
  fwrite(merged,file=signature_file,sep="\t")
}


# Assess co-localization of MPN and HSC signature -------------------------
library(BuenColors)
library(ggrastr)
library(ggridges)

signature_file <- "mpn_hsc_signatures file"
merged <- fread(signature_file)
cor.test(merged$hsc_signature,merged$mpn_signature,method="spearman")

cluster_sizes <- merged %>% group_by(cluster) %>% summarise(n = n())

longer <- pivot_longer(merged,cols=contains("signature"),names_to = "geneset",values_to = "score") %>%
  mutate(idx = paste(cluster,geneset,sep="-"))

# Check for significance of clusters
pvals <- lapply(unique(merged$cluster),function(clus){
  print(clus)
  subset <- merged %>% filter(cluster == clus)
  control <- merged %>% filter(cluster != clus)
  
  res1 <- wilcox.test(subset$hsc_signature,control$hsc_signature,alternative = "greater")
  res2 <- wilcox.test(subset$mpn_signature,control$mpn_signature,alternative = "greater")
  
  output <- data.frame(cluster = clus,
                       type = "HSC",
                       pval = res1$p.value)
  output <- rbind(output,data.frame(cluster = clus,
                                    type = "MPN",
                                    pval = res2$p.value))
  return(output)
}) %>% bind_rows()
pvals$fdr <- p.adjust(pvals$pval,method = "fdr")
pvals %>% filter(fdr < 0.001) %>% group_by(cluster) %>% 
  filter(n()  == 2) %>% .$cluster %>% unique()

pvals <- pvals %>% group_by(cluster) %>% 
  mutate(stars = case_when(max(fdr) < 0.001 ~ "**",
                   min(fdr) < 0.001 & max(fdr) > 0.001 ~ "*",
                   TRUE ~ "")) %>% ungroup()
cluster_sizes <- merge(cluster_sizes,unique(pvals[,c("cluster","stars")]),by="cluster") 

# Arrange in terms of decreasing signature cluster
data <- longer %>% filter(geneset=="mpn_signature") %>% group_by(cluster) %>% summarise(idx = mean(score)) %>% arrange(desc(idx))
longer$cluster <- factor(longer$cluster,levels = data$cluster)
cluster_sizes$cluster <- factor(cluster_sizes$cluster,levels = data$cluster)

# Set max score to 4
longer$score_capped <- MinMax(longer$score,min=min(longer$score),max=4)

dodge <- position_dodge(width = 0.5)
p1 <- ggplot(longer,aes(x=cluster))+
  # geom_violin(aes(group=idx,y=score_capped,fill=geneset)) +
  geom_boxplot_jitter(aes(group=idx,y=score_capped,fill=geneset),
                      outlier.size=0.05, outlier.jitter.width = 0.08, outlier.alpha=0.5,raster=T, raster.dpi = 400) +
  geom_text(data=cluster_sizes,aes(label=paste(n, stars)),y=4.2,size=1.5,angle = 90)+
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(2,3)]) +
  pretty_plot(fontsize=8) + L_border() +
  labs(x="Cluster",y="Score") + 
  theme(legend.position="none")
ggsave2(p1,file="../output/target_genes/scRNAseq/ICA_clusters_MPN_HSC_scores.pdf", height = 2,width =4.5)


# Group by HSC clusters 12,28,17 vs. all others
longer <- longer %>% mutate(HSC_binary = ifelse(cluster %in% c(12,28,17),"HSC","all others"))
longer$HSC_binary <- factor(longer$HSC_binary,levels = c("HSC","all others"))
p2 <- ggplot(longer,aes(x=HSC_binary))+
  # geom_boxplot(aes(group=paste0(HSC_binary,geneset))) +
  geom_boxplot_jitter(aes(group=paste0(HSC_binary,geneset),y=score_capped,fill=geneset),
                      outlier.size=0.05, outlier.jitter.width = 0.08, outlier.alpha=0.5,
                      raster=T, raster.dpi = 400) +
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(2,3)]) +
  pretty_plot(fontsize=8) + L_border() +
  labs(x="",y="Score") +
  theme(legend.position = "none")

ggsave2(p2,file="../output/target_genes/scRNAseq/ICA_HSC_vs_allothers.pdf", height = 2,width =1.5)

