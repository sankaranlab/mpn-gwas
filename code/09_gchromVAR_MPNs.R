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
set.seed(1026)
"%ni%" <- Negate("%in%")

trait <- "MPN_arraycovar_meta_finngen_r4"
towrite <- TRUE
# Import and run
# ATAC
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
                       list.files("../data/abf_finemap/for_gchromVAR/all", full.names = TRUE, pattern = "^MPN_array"),
                       colidx = 4)


# Run g-chromVAR
bg <- getBackgroundPeaks(SE,niterations=200)
dev <- computeWeightedDeviations(SE, mpns, background_peaks = bg)

# Reformat results
zscoreWeighted <- melt(t(assays(dev)[["z"]]))
zscoreWeighted[,2] <- gsub("_abf_cojo_PP0.001", "", zscoreWeighted[,2])
colnames(zscoreWeighted) <- c("Celltype","Trait","zscore")

zscoreWeighted$logp <- -log10(pnorm(zscoreWeighted$zscore, lower.tail = FALSE))
outdf <- zscoreWeighted %>% filter(!grepl('duplicate', Trait)) %>% arrange(desc(logp))
outdf$Celltype <- factor(outdf$Celltype,levels=outdf$Celltype)

# Barplot of -log10(p-values)
p1 <- ggplot(outdf,aes(x = Celltype, y = logp)) +
  geom_bar(width = 1, aes(fill = Celltype), colour="black",
           stat = "identity", position = position_dodge(width=1))+
  pretty_plot(fontsize = 8) + L_border() +
  scale_fill_manual(values = ejc_color_maps) +
  labs(x = "", y = "g-chromVAR Enrichment (-log10 p)", fill = "") +
  scale_y_continuous(expand=c(0.05,0))+
  theme(legend.text=element_text(size=8),
        legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1))
p1 

if (towrite){
  cowplot::ggsave2(p1, filename = paste0("../output/gchromVAR/gchromVAR_abf_",trait,"_logp.pdf"), width = 3, height = 3)
}

# Export table
if (towrite){
  fwrite(outdf, paste0("../output/gchromVAR/gchromVAR_abf_",trait,"_zscores.txt"), sep = "\t")
}

# z scores
p1 <- ggplot(outdf,aes(x = Celltype, y = zscore)) +
  geom_bar(width = 1, aes(fill = Celltype), colour="black",
           stat = "identity", position = position_dodge(width=1))+
  pretty_plot(fontsize = 10) + L_border() + 
  scale_fill_manual(values = ejc_color_maps) +
  labs(x = "", y = "g-chromVAR Enrichment (z-score)", fill = "") + 
  scale_y_continuous(expand=c(0.05,0))+
  theme(legend.text=element_text(size=8),
        legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1))

p1
if (towrite){
  cowplot::ggsave(p1, filename = paste0("../output/gchromVAR/gchromVAR_abf_",trait,"_zscores.pdf"), width = 4.5, height = 3)
}