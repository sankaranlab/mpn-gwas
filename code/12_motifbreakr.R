library(data.table)
library(tidyverse)
library(GenomicRanges)
library(BuenColors)
library(matrixStats)
library(SummarizedExperiment)
library(Matrix)
library(preprocessCore)
library(qvalue)
library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg19)
register(MulticoreParam(2))

# Read in ABF CS
CS.abf.df <- fread("../data/abf_finemap/MPN_arraycovar_meta_finngen_r4_abf_cojo_PP0.001.bed") 
names(CS.abf.df) <- c("seqnames","start","end","PP","region","var","rsid")
CS.abf.df$ref <- str_split_fixed(CS.abf.df$var,"_",3)[,2]
CS.abf.df$alt <- str_split_fixed(CS.abf.df$var,"_",3)[,3]

# Exclude indels
ref_exclude <- sapply(CS.abf.df$ref,function(y){nchar(as.character(y))}) >1 
alt_exclude <- sapply(CS.abf.df$alt,function(y){nchar(as.character(y))}) >1 
CS.abf.df_noindels <- CS.abf.df[!Reduce("|", list(ref_exclude,alt_exclude)),]

# Format to bed file input
CS.abf.df_noindels %>% distinct(var,.keep_all = T) %>%
  mutate(start = start-1) %>% mutate(end= start+1) %>%
  mutate(score = 0, strand = "+",var =paste0("chr",gsub("_",":",var))) %>%
  dplyr::select(seqnames,start,end,var,score,strand) -> CS.abf.bed

write.table(CS.abf.bed,file="../data/abf_finemap/r4_MPN_CML_abf_cojo_PP0.001.motifbreakr.bed",
            quote = F, sep = "\t", col.names = F, row.names = F)


#import the BED file
snps.mb.frombed <- snps.from.file(file = "../data/abf_finemap/r4_MPN_CML_abf_cojo_PP0.001.motifbreakr.bed",
                                  search.genome = BSgenome.Hsapiens.UCSC.hg19,
                                  format = "bed")

data(motifbreakR_motif)
data("hocomoco")

mpn_mbreaker <- motifbreakR(snpList = snps.mb.frombed, filterp = TRUE,
                            pwmList = hocomoco,
                            threshold = 1e-3,
                            method = "ic",
                            bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                            BPPARAM = BiocParallel::bpparam())

mpn_mbreaker <- calculatePvalue(mpn_mbreaker)
mpn_mbreaker$var <- names(mpn_mbreaker)

# Save motifbreakr results
saveRDS(mpn_mbreaker, file="../output/transcription_factors/r4_mpn_motifbreakr.rds")


# Reformat mbreaker results -----------------------------------------------
mpn_mbreaker <- readRDS("../output/transcription_factors/r4_mpn_motifbreakr.rds")

# Plot anecdotes
pdf(file=paste0("../output/transcription_factors/motif_plot_rs524137.pdf"), width = 3, height = 4)
par(cex.main=0.8,mar=c(1,1,1,1))
plotMB(mpn_mbreaker, rsid = "chr9:135879542:C:T", effect="strong")
dev.off()

# Convert ID back to UKID format
mpn_mbreaker$var <- sapply(mpn_mbreaker$var, function(var){
  temp <- str_split_fixed(gsub("chr","",var),":",4)
  paste0(temp[1],":",temp[2],"_",temp[3],"_",temp[4])
})
names(mpn_mbreaker) <- NULL

# Merge with finemap PPs
mpn_mbreaker.df <- merge(unique(CS.abf.df[,c("var","PP","rsid")]),unique(as.data.frame(mpn_mbreaker)),by="var",
                         allow.cartesian=TRUE) 

write.table(mpn_mbreaker.df,file="../output/transcription_factors/r4_mpn_motifbreakr_PPmerged.txt",
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE) 
