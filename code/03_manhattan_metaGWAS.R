library(qqman)
library(tidyverse)
library(data.table)
library(BuenColors)
library(GenomicRanges)
library(regioneR)
library(ggrastr)
library(ggman)
library(cowplot)
source("03a_manhattan_function.R")
"%ni%" <- Negate("%in%")

# Load meta sumstats
sumstats <- "file path to MPN summary statistics"
filtered <- fread(sumstats)
filtered <- filtered %>% mutate(start=POS,end=POS)
filtered.gr <- filtered %>% GRanges()

if (FALSE){
  # Find sentinels
  window_width <- 500000
  sentinels <- NULL
  for (chr in 1:23){
    subset <- filtered %>% filter(CHR == chr)
    if (min(subset$pvalue) > 1e-6){
      next 
    }
    
    allsugg <- subset %>% filter(pvalue < 1e-6)
    hit <- subset[subset$pvalue == min(subset$pvalue),]
    sentinels <- rbind(sentinels,hit)
    
    blacklist <- allsugg %>% filter(POS > hit$POS -window_width, POS < hit$POS + window_width)
    while (nrow(blacklist) < nrow(allsugg)){
      new_subset <- subset %>% filter(RSID %ni% blacklist$RSID)
      hit <- new_subset[new_subset$pvalue == min(new_subset$pvalue),]
      
      # Add on new sentinel to list and increase blacklist 
      sentinels <- rbind(sentinels,hit)
      blacklist <- rbind(blacklist,allsugg %>% filter(POS > hit$POS -window_width, POS < hit$POS + window_width))
    }
  }
  sentinels
}

# Use COJO sentinels to define loci
cojo_sentinels <- fread("../data/abf_finemap/MPN_arraycovar_meta_finngen_r4_abf_cojo_PP0.001_annotated.bed")%>%
  filter(sentinel=="yes") %>%  .$rsid
sentinels <- filtered %>% filter(RSID %in% cojo_sentinels)
  
window_width <- 500000
genome_wide <- sentinels %>% filter(pvalue < 5e-8) %>% GRanges()
suggestive <- sentinels %>% filter(pvalue > 5e-8, pvalue < 1e-6) %>% GRanges()

genome_wide.window <- extendRegions(genome_wide,extend.start = window_width,extend.end=window_width)
suggestive.window <- extendRegions(suggestive,extend.start = window_width,extend.end=window_width)

idx <- findOverlaps(genome_wide.window,filtered.gr)
sigsnps <- filtered[idx@to,]$RSID %>% unique()

idx <- findOverlaps(suggestive.window,filtered.gr)
suggsnps <- filtered[idx@to,]$RSID %>% unique()

tohighlight <- filtered %>% filter(RSID %in% sigsnps) %>% mutate(group="genome-wide")
tohighlight <- rbind(tohighlight, filtered %>% filter(RSID %in% suggsnps) %>% mutate(group="suggestive"))

filtered$newCHR <- paste0("chr",filtered$CHR)
p1 <- ggman(filtered %>% filter(pvalue<0.01), snp = "RSID", bp = "POS", chrom = "newCHR", pvalue = "pvalue",
            sigLine = 7.30103,relative.positions = TRUE,pointSize = 0.8) + 
  pretty_plot(fontsize = 12) + 
  scale_color_manual(values=c("gray", "dimgray")) +
  ggtitle("") +
  scale_y_continuous(breaks=seq(0,120,10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p2 <- p1 +  geom_point_rast(data = p1[[1]] %>% filter(RSID %in% suggsnps),colour= "darkolivegreen3",size=0.8) +
  geom_point_rast(data = p1[[1]] %>% filter(RSID %in% sigsnps),colour= "blue",size=0.8) +
  geom_point_rast(data = p1[[1]] %>% filter(RSID %in% cojo_sentinels),colour= "black",fill="firebrick",pch=21,size=1.5) 

cowplot::ggsave2(p2, file="../output/gwas_plots/meta.finngen.r4.manhattan.pdf",
                width=9.5,height=5.5)


# QQ plot -------------------------------------------------------------
# Lambda GC
lambda.gc(filtered$pvalue)

# Q-Q plot
filtered$pvalue <- as.numeric(as.character(filtered$pvalue))
filtered <- filtered[complete.cases(filtered),]

pvector <- filtered$pvalue
toplot <- data.frame(observed = round(-log10(sort(pvector,decreasing=FALSE)),3),
                     expected = round(-log10( ppoints(length(pvector) )),3))
toplot.filtered <- unique(toplot)

p1 <- ggplot(toplot.filtered,aes(x=expected,y=observed))+
  geom_point_rast() +
  geom_abline(intercept=0,slope=1)+
  labs(x="Expected -log10(p)",y="Observed -log10(p)") + 
  pretty_plot(fontsize = 8) + L_border()
cowplot::ggsave2(p1,filename="../output/gwas_plots/meta.finngen.r4.qq.pdf",height=4,width=4)

