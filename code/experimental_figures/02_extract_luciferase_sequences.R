library(tidyverse)
library(data.table)
library(DNAshapeR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(preprocessCore)
library(matrixStats)
library(seqinr)

# Variants of interest
# PODXL: 7:130742066_A_G	rs7803075
# MECOM: 3:168822748_G_T rs13327022 
# FOXO1: 13:41204015_T_C	rs7323267
# TERT: 5:1138335_T_G	rs4131149
# RUNX1: 21:36351891_T_C	rs2834712

vars <- c("rs7803075","rs13327022","rs7323267","rs4131149","rs2834712")
  
# Import CS 
CS.df <- fread("../data/abf_finemap/MPN_CML_abf_cojo_95CS.bed")

# 500bp peaks/counts
peaksdf <- fread("../data/atac/29August2017_EJCsamples_allReads_500bp.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts.df <-  data.matrix(fread("../data/atac/29August2017_EJCsamples_allReads_500bp.counts.txt"))

# Log2 cpm normalize
cpm <- round(sweep(counts.df, 2, colSums(counts.df), FUN="/") * 1000000, 1)
log2cpm <- log2(cpm+1)
# Min / max scale
log2cpm.minmax <- log2cpm / rowMax(log2cpm)
# Select for myeloid populations
myeloid.counts <- as.data.frame(log2cpm) %>% 
  dplyr::select("HSC", "MPP", "CMP", "MEP", "Ery","Mega") %>% as.matrix()

  
luc <- CS.df %>% filter(rsid %in% vars) %>% mutate(pos = start)
luc$alt <- str_split_fixed(luc$var,"_",3)[,3]

# Overlap
idx <- findOverlaps(peaks,GRanges(luc))
fm.counts <- myeloid.counts[idx@from,]

# Construct mega data.frame with peaks, cpm, and finemap info
mega_df <- data.frame(
  luc[idx@to,c("var","rsid","region","PP","pvalue","region_rank","pos","alt")],
  peaks[idx@from],
  myeloid.counts[idx@from,]
) %>% dplyr::select(-strand,-width) %>% unique()

mega_df$dist <- mega_df$pos - mega_df$start +1


# Create own GRanges file
# 7:130741947-130742264
# 3:168822674-168822855
# 13:41203850-41204167
# 5:1138148-1138578
# 21:36351658-36352066

luc_seqs <- data.frame(
  var=c("7:130742066_A_G","3:168822748_G_T","13:41204015_T_C","5:1138335_T_G","21:36351891_T_C"),
  seqnames = paste0("chr",c(7,3,13,5,21)),
  start=c(130741947,168822674,41203850,1138148,36351658),
  end=c(130742264,168822855,41204167,1138578,36352066)
)

# Plot GFI1B fasta
luc_seqs <- data.frame(
  var=c("9:135870130_C_G"),
  seqnames = paste0("chr",c(9)),
  start=c(135870130-50),
  end=c(135870130+50)
)
luc_seqs$pos <- as.integer(gsub("_.*","",str_split_fixed(luc_seqs$var,":",2)[,2]))
luc_seqs$alt <- str_split_fixed(luc_seqs$var,"_",3)[,3]
luc_seqs$dist <- luc_seqs$pos - luc_seqs$start +1

GRanges(luc_seqs)

alt_seq <- ref_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19,GRanges(luc_seqs),as.character=TRUE) 
write.fasta(as.list(ref_seq),paste0(luc_seqs$var,"_ref"),file.out="../../output/luciferase_sequences/luc_variant_sequences.fa")

# Modify one variant in alternate sequence
sapply(seq(1,length(alt_seq)),function(y){
  dist <- luc_seqs$dist[y]
  alt <-  luc_seqs$alt[y]
  substring(alt_seq[y], dist, dist) <- alt
  return(alt_seq[y])
}) -> alt_seq_changed

write.fasta(as.list(alt_seq_changed),paste0(luc_seqs$var,"_alt"),open="a",
            file.out="../../output/luciferase_sequences/luc_variant_sequences.fa")
