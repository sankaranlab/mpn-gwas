library(chromVAR)
library(gchromVAR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(tidyverse)
library(SummarizedExperiment)
library(Matrix)
library(BuenColors)
library(cowplot)
library(diffloop)
library(preprocessCore)
library(stringr)
"%ni%" <- Negate("%in%")

set.seed(1026)

# Import and run run run
# # 250bp
peaksdf <- fread("../data/atac/26August2017_EJCsamples_allReads_250bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <-  data.matrix(fread("../data/atac/26August2017_EJCsamples_allReads_250bp.counts.txt"))
cpm <- round(sweep(counts, 2, colSums(counts), FUN="/") * 1000000, 1)

# Create objects for g-chromVAR
SE <- SummarizedExperiment(assays = list(counts = counts),
                           rowData = peaks, 
                           colData = DataFrame(names = colnames(counts)))
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)

# Function to take the sum of the absolute values of a number
# Needs at least 2 bed files to compare with each other
mpns <- importBedScore(rowRanges(SE), 
                       list.files("../data/abf_finemap/for_gchromVAR/all", full.names = TRUE, pattern = ".bed$"),
                       colidx = 4)

# Analyze where weights are concentrated
if (FALSE){
  keepPeaks <- assays(mpns)[["weights"]][,1] > 0.01
  mega_df <- data.frame(
    rowRanges(SE)[keepPeaks],
    data.matrix(cpm[keepPeaks,]),
    sumPP=data.matrix(assays(mpns)[["weights"]][keepPeaks,1]))
  mega_df %>% arrange(desc(sumPP)) %>% dplyr::select(seqnames,start,end,HSC,CLP,sumPP)
}

# Run g-chromVAR
bg <- getBackgroundPeaks(SE,niterations=50)
dev <- computeWeightedDeviations(SE, mpns, background_peaks = bg)
outdf <- melt(t(assays(dev)[["z"]]))
outdf$pval <- pnorm(outdf$value, lower.tail = FALSE)
outdf$logp <- -log10(outdf$pval)
mdf <- outdf %>% filter(!grepl('duplicate', Var2)) %>% arrange(desc(logp))
mdf$qvalue <- qvalue::qvalue(mdf$pval,lambda=0)$qvalues
mdf$Var1 <- factor(mdf$Var1,levels=mdf$Var1)
mdf

# Barplot of -log10(p-values)
p1 <- ggplot(mdf,aes(x = Var1, y = logp)) +
  geom_bar(width = 1, aes(fill = Var1), colour="black",
           stat = "identity", position = position_dodge(width=1))+
  pretty_plot(fontsize = 8) + L_border() +
  scale_fill_manual(values = ejc_color_maps) +
  labs(x = "", y = "g-chromVAR Enrichment (-log10 p)", fill = "") +
  scale_y_continuous(expand=c(0.05,0))+
  theme(legend.text=element_text(size=8),
        legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1))
p1 

if (FALSE){
  cowplot::ggsave2(p1, filename = "../output/gchromVAR/gchromVAR_abf_MPN_noncoding_only_logp.pdf", width = 3, height = 3)
}

# Make tables and export
widemat <- t(assays(dev)[["z"]])
colnames(widemat) <- gsub(".PP001","",colnames(widemat))
zscoreCHROMVAR <- melt(widemat)
zscoreCHROMVAR <- zscoreCHROMVAR %>% filter(!grepl('duplicate', Var2)) %>% arrange(desc(abs(value)))
write.table(zscoreCHROMVAR, "../output/gchromVAR/gchromVAR_abf_MPN_zscores.txt", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# z scores
p1 <- ggplot(mdf,aes(x = Var1, y = value)) +
  geom_bar(width = 1, aes(fill = Var1), colour="black",
           stat = "identity", position = position_dodge(width=1))+
  pretty_plot(fontsize = 10) + L_border() + 
  scale_fill_manual(values = ejc_color_maps) +
  labs(x = "", y = "g-chromVAR Enrichment (z-score)", fill = "") + 
  scale_y_continuous(expand=c(0.05,0))+
  theme(legend.text=element_text(size=8),
        legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1))

p1
if (FALSE){
  cowplot::ggsave2(p1, filename = "../output/gchromVAR/gchromVAR_abf_MPN_zscores.pdf", width = 4.5, height = 3)
}


