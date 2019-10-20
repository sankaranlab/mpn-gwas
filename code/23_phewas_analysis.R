library(data.table)
library(tidyverse)
library(qvalue)
library(Matrix)
library(ComplexHeatmap)
library(SummarizedExperiment)
library(BuenColors)
library(annotables)
library(ggrastr)
library(ggrepel)
"%ni%" <- Negate("%in%")

# Write list of variants for pheWAS input ---------------------------------
if (FALSE){
  # Add leading 0 to chromosomes under 10
  CS.df <- fread("../data/abf_finemap/MPN_CML_abf_cojo_PP0.001_annotated.bed") %>% 
    mutate(chr = as.integer(gsub("chr","",seqnames))) %>% mutate(var = ifelse(chr < 10,paste0("0",var),var)) %>%
    dplyr::select(var) %>% unique()
  fwrite(CS.df,file="../output/phewas/MPN_PP001_phewas_input_leading_zeros.txt",sep="\t",col.names = F)
}

# SAIGE pheWAS ------------------------------------------------------------
# Read phewas results
icd <- fread("zcat < ../output/phewas/SAIGE.allphenos.MPN.txt.gz") %>%
  mutate(var = paste(CHROM,POS,REF,ALT,sep=":"), maf = ifelse(af < 0.5, af, 1-af),
         expected_case_minor_AC =  2 * maf * num_cases, logp = -log10(pval),
         var = paste0(CHROM,":",POS,"_",REF,"_",ALT)) %>%
  filter(expected_case_minor_AC > 50)

# Read and merge pheno codes 
phecodes <- readxl::read_xlsx("../output/phewas/saige-phenotype-information.xlsx") %>% 
  dplyr::rename(pheno = `Phenotype Description`,category=`Phenotype Category`) %>%
  filter(category %ni% c("pregnancy complications","injuries & poisonings","sense organs","symptoms"))
phecodes$PheCode <- as.double(phecodes$PheCode)
icd <- merge(icd,phecodes[,c("PheCode","pheno","category")],by.x="phecode",by.y="PheCode")

# Combine with sum stats 
CS.df <- fread("../data/abf_finemap/MPN_CML_abf_cojo_95CS.bed")
merged <- merge(icd,CS.df[,c("var","PP","region_rank")],by="var") %>%
  mutate(CHROM=paste0("chr",CHROM))

# Remove MPNs and polycythemia vera due to redundancy
merged <- merged %>% filter(pheno %ni% c("Myeloproliferative disease","Polycythemia vera","Polycythemia, secondary"))

# Exclude HLA region
hla_bounds <- c(28866528,33775446)
toexclude <- merged %>% filter(CHROM=="chr6", POS > hla_bounds[1],POS < hla_bounds[2]) %>% .$var
merged <- merged %>% filter(var %ni% toexclude)

# Add nearest gene
# Read in gene body annotations
gencode_gr <- readRDS("../data/annotations/gencode_filtered.rds")

# Find nearest gene
merged.gr <- makeGRangesFromDataFrame(merged,keep.extra.columns = T,seqnames.field="CHROM",
                                      start.field = "POS",end.field = "POS")
merged$nearest_gene <- gencode_gr[nearest(merged.gr,gencode_gr,ignore.strand=TRUE),]$gene_name

# Set bonferroni threshold
highPP <- merged %>% filter(PP > 0.10| region_rank ==1 )
num_variants <- highPP %>% .$var %>% unique %>% length() 
p_threshold <- 0.05 / (length(unique(merged$pheno)))

tolabel <- highPP %>% group_by(category) %>% 
  dplyr::slice(which.max(logp)) %>%
  mutate(tolabel=ifelse(logp > -log10(p_threshold),paste(pheno,nearest_gene,sep="-"),"")) %>%  
  arrange(category,desc(logp)) %>% as.data.frame() %>% dplyr::select(var,tolabel,pheno) 
toplot<- right_join(tolabel,highPP,by=c("var","pheno")) %>%
  mutate(tolabel = ifelse(is.na(tolabel),"",tolabel)) %>% 
  group_by(pheno) %>%mutate(mx = max(logp)) %>%
  arrange(category,desc(mx),desc(logp)) %>% ungroup() %>%   mutate(order = row_number())
toplot$pheno <- factor(toplot$pheno, levels=unique(toplot$pheno))

toplot$FDR <- qvalue::qvalue(toplot$pval)$qvalues

p1 <- ggplot(toplot,aes(x=pheno,y=-log10(FDR),label=tolabel))+
  geom_point_rast(aes(color=category,fill=category),alpha=1,size=0.5) +
  scale_color_manual(values=jdb_palette("lawhoops"))+
  geom_hline(yintercept = -log10(0.01), linetype = 2) +
  labs(y="pheWAS -log10(FDR)",x="Phenotype") +
  geom_text_repel(angle = 0,size=2) +  
  pretty_plot(fontsize=8) + L_border() +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(),
        strip.background = element_blank(), strip.text = element_blank(),
        legend.position ="bottom") +
  scale_y_continuous(expand = c(0.05, 0))
p1
cowplot::ggsave2(p1,file="../output/phewas/phewas_dotplot_categories.pdf",width=4,height=3)


# Write table
phewas_table <- merged %>%filter(pval < p_threshold) %>% filter(PP > 0.10| region_rank ==1 )%>%
  dplyr::select(var,ID,pheno,num_cases,maf,beta,sebeta,pval,PP,region_rank,nearest_gene) %>% 
  arrange(pval)
fwrite(phewas_table,file="../output/phewas/phewas_sig_phenos.tsv",sep = "\t")

