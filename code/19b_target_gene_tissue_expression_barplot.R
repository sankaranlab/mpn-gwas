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

target_genes <- as.character(read.table("../output/target_genes/r4_target_genes_scored_noABC_ties.txt")$V1)

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
# Rank sum test
log2cpm.df <- data.frame(log2cpm) %>% rownames_to_column(var = "gene") %>%
  mutate(target_gene = gene %in% target_genes)
permuted <- sapply(1:10000, function(i) 
  sum(1:nrow(log2cpm.df) * sample(log2cpm.df$target_gene, length(log2cpm.df$target_gene))))

celltypes <- colnames(log2cpm) %>% gsub("-",".",.)
ranksums <- sapply(celltypes,function(ct){
    sum(1:nrow(log2cpm.df)*log2cpm.df[order(log2cpm.df[,ct], decreasing = TRUE), "target_gene"])
  })

pvals <-  as.numeric(2*pnorm((mean(permuted) - ranksums)/sd(permuted), lower.tail = FALSE)) %>%
  data.frame() %>% rownames_to_column() %>% setNames(.,c("celltype","pval")) %>% 
  mutate(celltype = colnames(log2cpm),FDR = p.adjust(pval,method = "fdr"),logp = -log10(pval)) 

pvals
pvals$celltype <- factor(pvals$celltype,levels = colnames(tg))
p2 <- ggplot(pvals,aes(x=celltype,y=-log10(pval)))+
  geom_bar(width = 1, aes(fill = celltype), colour="black",stat = "identity") +
  scale_fill_manual(values = ejc_color_maps) +
  geom_hline(yintercept = -log10(0.05/ncol(log2cpm)),linetype="dashed") +
  pretty_plot(fontsize=8) + L_border() +
  theme(legend.position = "none",axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +scale_y_continuous(expand = c(0, 0),limits=c(0,5.5)) 
p2

if (FALSE){
  cowplot::ggsave2(p2,file="../output/target_genes/r4_target_gene_scored_noABC_min2_ranksum_ct_enrichment.barplot.pdf",width=2.99,height=0.8)
  fwrite(pvals,file="../output/target_genes/r4_target_gene_scored_noABC_min2_ranksum_ct_enrichment.tsv",sep="\t")
}

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

box <- ggplot(toplot, aes(x=variable, y=value,fill=variable, color=variable)) + 
  geom_boxplot(varwidth = F,alpha = 0.6) + 
  scale_fill_manual(values = ejc_color_maps) + 
  scale_color_manual(values = ejc_color_maps) + 
  geom_point(data = combined, aes(x=variable, y=pc_mean),shape = 23,color="black",fill="black", size = 1.5) +
  pretty_plot(fontsize=8) + L_border() + 
  labs(x="",y="Target gene expression (log2(cpm))") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) 
box

p3 <- ggplot(combined, aes(x=variable, y=normalized,fill=variable)) + 
  geom_bar(position=position_identity(), stat="identity", color = "black", width = 0.7) +
  scale_fill_manual(values = ejc_color_maps) +
  labs(x="",y="Normalized target gene expression") +
  pretty_plot(fontsize=8) + L_border() +
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) 
p3

if (FALSE){
  cowplot::ggsave2(p2, file="../output/target_genes/r4_scored_noABC_min2_target_genes_expression_barplot_allheme.pdf",
                   height = 2.5,width = 3)
  cowplot::ggsave2(box, file="../output/target_genes/r4_scored_noABC_min2_target_genes_expression_boxplot_allheme.pdf",
                   height = 2.5,width = 3)
  cowplot::ggsave2(p3, file="../output/target_genes/r4_scored_noABC_min2_target_genes_expression_barplot_allheme_normalized.pdf",
                   height = 2.5,width = 3)
}
