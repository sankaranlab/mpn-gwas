library(data.table)
library(tidyverse)
library(BuenColors)
library(scales)
library(preprocessCore)
library(annotables)
library(qvalue)
"%ni%" <- Negate("%in%")
set.seed(1026)

eryth_color_maps <- c("P1" = "#3b82ae", "P2" = "#547294", "P3" = "#6d617a", "P4" = "#865160", "P5" = "#9f4046", "P6" = "#b8302c", "P7" = "#d11f12", "P8" = "#de1705")
all_color_maps <- c(ejc_color_maps[c("HSC", "MPP", "CMP")], "MEP" = "#FF81AF" ,eryth_color_maps)

target_genes <- read.table("../output/target_genes/target_genes_combined.PP10.txt")$V1

# All Hematopoiesis -------------------------------------------------------
# Read in Corces RNA data
rna <- data.frame(fread(paste0("../data/rna/16populations_RNAcounts.txt")))
genes <- rna[,1]
counts <- data.matrix(rna[,-1])
cpm <- sweep(counts, 2, Matrix::colSums(counts), FUN="/") * 1000000
colnames(cpm) <- c("HSC", "MPP", "LMPP", "CLP", "GMP-A", "GMP-B", "GMP-C", "CMP", "MEP", "NK", "CD4", "CD8", "B", "pDC", "Mono", "Ery")
rownames(cpm) <- genes

# Filter genes
keepGenes <- which(rowSums(counts) >= 0)
length(keepGenes)
cpm <- cpm[keepGenes, ]
log2cpm <- log2(cpm+1)

# Filter for protein coding genes
grch38.pc <- grch38 %>% filter(biotype == "protein_coding")
log2cpm <- log2cpm[rownames(log2cpm) %in% grch38.pc$symbol,]
log2cpm <- log2cpm[rowSums(log2cpm) > 0,]

# Calculate cross-tissue enrichments
quantile_normalize <- FALSE
if (quantile_normalize){
  qnorm_log2cpm <- normalize.quantiles(log2cpm)
} else{
  qnorm_log2cpm <- log2cpm
}
rownames(qnorm_log2cpm) <- rownames(log2cpm); colnames(qnorm_log2cpm) <- colnames(log2cpm)

tg <- qnorm_log2cpm[rownames(qnorm_log2cpm) %in% target_genes,]

# Rank sum test of target genes in each cell type
apply(tg,1,rank)%>% rowSums()

# Z score within each cell type, meta
ct_zcores <- sapply(target_genes,function(y){
  subset <- log2cpm[rownames(log2cpm)==y,]
  sds <- apply(log2cpm,2,sd)
  zscore <- (subset - colMeans(log2cpm)) / sds
  return(zscore)
})
meta_ct_zcore <- rowSums(ct_zcores) / sqrt(ncol(ct_zcores))
meta_ct_zcore
enrich_pvals <- -log10(pnorm(meta_ct_zcore,lower.tail = FALSE)) %>% data.frame() %>% rownames_to_column() %>% setNames(.,c("celltype","logp"))

ggplot(enrich_pvals,aes(x=celltype,y=logp))+
  geom_bar(width = 1, aes(fill = celltype), colour="black",stat = "identity") +
  scale_fill_manual(values = ejc_color_maps) +
  geom_hline(yintercept = -log10(0.05/ncol(log2cpm))) + 
  pretty_plot(fontsize=8) + L_border()


# Rank sum test
log2cpm.df <- data.frame(log2cpm) %>% rownames_to_column(var = "gene") %>%
  mutate(target_gene = gene %in% target_genes)
permuted <- sapply(1:10000, function(i) 
  sum(1:nrow(log2cpm.df) * sample(log2cpm.df$target_gene, length(log2cpm.df$target_gene))))

celltypes <- colnames(log2cpm) %>% gsub("-",".",.)
ranksums <- sapply(celltypes,function(ct){
    sum(1:nrow(log2cpm.df)*log2cpm.df[order(log2cpm.df[,ct], decreasing = TRUE), "target_gene"])
  })

pvals <- pnorm((mean(permuted) - ranksums)/sd(permuted), lower.tail = FALSE)%>% data.frame() %>% rownames_to_column() %>% setNames(.,c("celltype","pval")) %>% 
  mutate(celltype = colnames(log2cpm),FDR = p.adjust(pval,method = "fdr"),logp = -log10(pval)) 

pvals
pvals$celltype <- factor(pvals$celltype,levels = colnames(tg))
p1 <- ggplot(pvals,aes(x=celltype,y=-log10(pval)))+
  geom_bar(width = 1, aes(fill = celltype), colour="black",stat = "identity") +
  scale_fill_manual(values = ejc_color_maps) +
  # geom_hline(yintercept = -log10(0.05/ncol(log2cpm)),linetype="dashed") + 
  pretty_plot(fontsize=8) + L_border() +
  theme(legend.position = "none",axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +scale_y_continuous(expand = c(0, 0))
p1

if (FALSE){
  cowplot::ggsave2(p1,file="../output/target_genes/target_gene_ranksum_ct_enrichment.barplot.pdf",width=2.99,height=0.8)
  fwrite(pvals,file="../output/target_genes/target_gene_ranksum_ct_enrichment.tsv",sep="\t")
}


zscores <- apply(tg, 1, scale)
rownames(zscores) <- colnames(tg)
meta_z <- rowSums(zscores) / sqrt(ncol(zscores))
meta_z

enrich_pvals <- pnorm(meta_z,lower.tail = FALSE) %>% data.frame() %>% rownames_to_column() %>% setNames(.,c("celltype","pval"))
enrich_pvals

fwrite(enrich_pvals,file="../output/target_genes/target_gene_stouffer_z_ct_enrichment.tsv",sep="\t")

t.test(as.numeric(tg[,c("HSC","MPP")]),
            tg[,c("NK","CD4","CD8","B","pDC","Mono","Ery")])


# Expression barplots -----------------------------------------------------
# Plot target gene cpm barplot
log2cpm.df <- as.data.frame(qnorm_log2cpm)
log2cpm.df$gene <- rownames(log2cpm.df) 
log2cpm_long <- melt(log2cpm.df)

toplot <- log2cpm_long %>% filter(gene %in% target_genes)
tg_cpm_summarized <- toplot %>%
  group_by(variable) %>% summarize(tg_mean = mean(value),tg_se = sd(value)/sqrt(n())) 
all_pc <- log2cpm_long %>%  filter(gene %ni% target_genes) %>%
  group_by(variable) %>% summarize(pc_mean = mean(value),pc_se = sd(value)/sqrt(n()))

# Normalize by all other protein coding genes
combined <- merge(tg_cpm_summarized,all_pc,by="variable")
combined <- combined %>% mutate(normalized = tg_mean/pc_mean)
combined

p1 <- ggplot(data=toplot,aes(x=variable,y=value,group=gene)) + 
  geom_line(alpha=0.5,color="dodgerblue")+
  pretty_plot() + L_border()

p2 <- ggplot(combined, aes(x=variable, y=tg_mean,fill=variable)) + 
  geom_bar(position=position_identity(), stat="identity", color = "black", width = 0.7) +
  geom_point(aes(y=pc_mean)) +
  scale_fill_manual(values = ejc_color_maps) + 
  geom_errorbar(aes(ymin=tg_mean-tg_se, ymax=tg_mean+tg_se),
                width=0.3, position=position_dodge(.9)) +
  labs(x="",y="Target gene expression (log2(cpm))") +
  pretty_plot(fontsize=8) + L_border() +
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) 
p2

cowplot::ggsave2(p2, file="../output/target_genes/target_genes_expression_barplot_allheme.pdf",
                height = 2.5,width = 3)

p3 <- ggplot(combined, aes(x=variable, y=normalized,fill=variable)) + 
  geom_bar(position=position_identity(), stat="identity", color = "black", width = 0.7) +
  scale_fill_manual(values = ejc_color_maps) +
  labs(x="",y="Normalized target gene expression") +
  pretty_plot(fontsize=8) + L_border() +
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) 
p3

cowplot::ggsave2(p3, file="../output/target_genes/target_genes_expression_barplot_allheme_normalized.pdf",
                height = 2.5,width = 3)


