library(tidyverse)
library(data.table)
library(BuenColors)
library(cowplot)
library(GenomicRanges)
library(diffloop)
"%ni%" <- Negate("%in%")

toplot <- T

# Read in CS and sentinels
CS.df <- fread("../data/abf_finemap/MPN_CML_abf_cojo_95CS.bed")
sentinels <- fread("../data/sentinels/MPN_CML.meta.cojo.suggestive_loci.withsentinels.tsv")

# Median genetic width of CS
CS.df %>% group_by(region) %>% summarise(ranges=max(start)-min(start)) %>% .$ranges %>% median()

# Correlate PP with p-value, grouped by region 
pval_cors <- lapply(seq(1,length(unique(CS.df$region))),function(y){
  filtered <- CS.df %>% filter(region == y)
  ggplot(filtered,aes(x=-log10(pvalue),y=PP)) +
    geom_point()+
    ylim(0,1) + 
    pretty_plot() + L_border()
})

# Size of CS
sentinels$CSbin <- cut(sentinels$CS, c(1,5,10,20,50,10000),include.lowest=T)
CS.df.sum <- sentinels %>%
  group_by(CSbin) %>%
  summarize(count=n()) 
p1 <- ggplot(CS.df.sum,aes(y=count,x=1)) + 
  geom_bar(stat="identity",aes(fill=CSbin),position = position_stack(reverse = TRUE)) + 
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 3) +
  coord_flip() +
  pretty_plot(fontsize = 8) + L_border()+
  scale_y_continuous(expand=c(0,0))+
  labs(x="",y="number of regions")+
  scale_fill_manual(values = c(jdb_palette("GrandBudapest2")[c(3,4)],jdb_palette("Zissou")[1:5])) +
  theme_void() +
  theme(legend.position="none") +
  guides(fill=guide_legend(title="CS size"))

if (toplot){
  cowplot::ggsave22(p1, file="../output/finemap/number_variants_each_CS.pdf", width=6.5,height=1)
}

# Correlate size of CS with MAF
p1 <- ggplot(sentinels,aes(x=freq,y=CS))+
  geom_point(color="blue",alpha=0.5)+
  pretty_plot(fontsize=8) + L_border() +
  labs(x="sentinel MAF",y="size of credible set")

p2 <- ggplot(sentinels,aes(x=freq,y=CS))+
  geom_point(color="blue",alpha=0.5)+
  scale_y_continuous(trans='log10') +
  pretty_plot(fontsize=8) + L_border() +
  labs(x="sentinel MAF",y="size of credible set")

if (toplot){
  cowplot::ggsave2(p1, file="../output/finemap/maf_to_CS_size.pdf", width=2,height=2)
  cowplot::ggsave2(p2, file="../output/finemap/maf_to_CS_size_log10.pdf", width=2,height=2)
}

# Bin the CS PPs for high threshold
CS.df.highPP <- CS.df[CS.df$PP> 0.25,]
CS.df.highPP$PPbin <- cut(CS.df.highPP$PP, c(0.250,0.5,0.75,0.90,0.95,0.99,1.00))
CS.df.sum <- CS.df.highPP %>%
  group_by(PPbin) %>%
  summarize(count=n()) 

p2 <- ggplot(CS.df.sum,aes(y=count,x=1)) + 
  geom_bar(stat="identity",aes(fill=PPbin),position = position_stack(reverse = T)) + 
  geom_text(aes(label = count), position = position_stack(vjust = 0.5,reverse=F), size = 3) +
  coord_flip()+
  pretty_plot(fontsize = 8) + L_border()+
  scale_y_continuous(expand=c(0,0))+
  labs(x="",y="number of regions")+
  scale_fill_manual(values = c(jdb_palette("GrandBudapest2")[c(3,4)],jdb_palette("Zissou")[1:5])) +
  theme_void() +
  theme(legend.position="none")s

# Top variant PP per region
CS.best <- CS.df %>% arrange(region,desc(PP)) %>% group_by(region) %>% top_n(1,PP)
CS.best$PPbin <- cut(CS.best$PP, c(0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.75, 1.0))
CS.df.best.sum <- CS.best %>%
  group_by(PPbin) %>%
  summarize(count=n())  %>% mutate(trait="MPN")
CS.df.best.sum$PPbin <- factor(CS.df.best.sum$PPbin, levels=levels(CS.df.best.sum$PPbin))

p4 <- ggplot(CS.df.best.sum,aes(y=count,x=1,fill=PPbin,label=count)) + 
  geom_col() + 
  geom_text(position = position_stack(vjust = 0.5), size = 3) +
  coord_flip() +
  pretty_plot(fontsize = 8) + L_border()+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values = c(jdb_palette("GrandBudapest2")[c(4)],jdb_palette("Zissou")[1:5]))+
  guides(fill = guide_legend(reverse = FALSE))+
  theme_void() +
  theme(legend.position="none") 
if (toplot){
  cowplot::ggsave2(p4, file="../output/finemap/top_variant_per_region.pdf", width=6.5,height=1)
}

# Genomic annotations
x01_coding_gr <- bedToGRanges("../data/annotations/Coding_UCSC.bed")
x02_promoter_gr <- bedToGRanges("../data/annotations/Promoter_UCSC.fixed.bed")
x03_utr_gr <- bedToGRanges("../data/annotations/UTR_3_UCSC.bed")
x04_atac_gr <- bedToGRanges("../data/atac/26August2017_EJCsamples_allReads_250bp.bed")
x05_intron_gr <- bedToGRanges("../data/annotations/Intron_UCSC.bed")

CS.gr <- CS.df %>% filter(PP> 0.01) %>% GRanges()
end(CS.gr) <- start(CS.gr)

# Do overlaps
ov_1 <- findOverlaps(CS.gr, x01_coding_gr)
ov_2 <- findOverlaps(CS.gr, x02_promoter_gr)
ov_3 <- findOverlaps(CS.gr, x03_utr_gr)
ov_4 <- findOverlaps(CS.gr, x04_atac_gr)
ov_5 <- findOverlaps(CS.gr, x05_intron_gr)

# Classify each accessibility peak
classAll <- ifelse(1:length(CS.gr) %in% queryHits(ov_1), "coding",
                   ifelse(1:length(CS.gr) %in% queryHits(ov_2), "promoter",
                          ifelse(1:length(CS.gr) %in% queryHits(ov_3), "utr",
                                 ifelse(1:length(CS.gr) %in% queryHits(ov_4), "accessible",
                                        ifelse(1:length(CS.gr) %in% queryHits(ov_5), "intron", "intergenic")))))

table(classAll)*100/sum(table(classAll))
order <- c("coding", "promoter", "accessible", "intron","intergenic")
totalDF <- data.frame(
  what = order,
  prop = as.numeric(table(classAll)[order]*100/sum(table(classAll)))
) 

totalDF$what <- factor(totalDF$what,levels=order)

p5 <- ggplot(totalDF, aes(x = 1, y = prop)) + 
  geom_bar(stat="identity",aes(fill=what),position = position_stack(reverse = TRUE)) + 
  geom_text(aes(label = round(prop,1)), position = position_stack(vjust = 0.5 ), size = 3) +
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(3,2,1,5,8)])+
  scale_y_continuous(expand=c(0,0))+
  coord_flip()+
  pretty_plot() + L_border() +
  theme_void() +
  theme(legend.position = "none",legend.title = element_blank())
  
if (toplot){
  cowplot::ggsave2(p5, filename = "../output/finemap/annotations_proportions.pdf", width = 6.5,height = 1)
}
