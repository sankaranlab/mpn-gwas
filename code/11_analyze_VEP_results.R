library(tidyverse)
library(data.table)
library(BuenColors)

# Read in CS
CS.df <- fread("../data/abf_finemap/MPN_CML_abf_cojo_PP0.001_annotated.bed") %>% 
  mutate(ref = str_split_fixed(var,"_",3)[,2],alt = str_split_fixed(var,"_",3)[,3]) %>%
  group_by(var) %>% dplyr::slice(which.max(PP)) %>% mutate(pos = start) %>% ungroup()

# Write table for VEP
if (FALSE){
  for_vep <- CS.df  %>%dplyr::select(seqnames,pos,ref,alt) %>% 
    mutate(seqnames = gsub("chr","",seqnames),filler = ".",filler2=".",filler3=".",filler4=".") %>%
    dplyr::select(seqnames,pos,filler,everything())%>% unique()
  
  write.table(for_vep,file="../output/VEP/gcta_cojo_PP001.vcf",quote = FALSE, sep = "\t", col.names = F, row.names = F)
}

# Read in VEP
vep_mostsevere <- fread("../output/VEP/VEP_MPN_PP001_most_severe_consequence.txt") %>% 
  dplyr::rename("most_severe_consequence"=Consequence) %>% dplyr::select(Location,most_severe_consequence)
vep_all <- fread("../output/VEP/VEP_MPN_PP001.txt") %>% 
  dplyr::select(Location,BIOTYPE,SYMBOL,EXON,SIFT,PolyPhen,CADD_PHRED) %>% arrange(Location,desc(SIFT)) %>% 
  distinct(Location,.keep_all = T) 
vep_combined <- merge(vep_all,vep_mostsevere,by="Location") %>% mutate(Location = gsub("-.*","",Location))

# Merge with VEP 
CS.df <- CS.df %>% mutate(Location = gsub("_.*","",var))
CS.df.vep <- as.data.frame(left_join(CS.df,vep_combined,by="Location")) %>% arrange(region,region_rank)
CS.df.vep.PP01 <- CS.df.vep %>% filter(PP>0.01)

# Filter for coding variants
table(CS.df.vep$most_severe_consequence)
coding_consequences <- c("missense_variant","synonymous_variant","frameshift_variant",
                         "splice_acceptor_variant","splice_donor_variant","splice_region_variant",
                         "inframe_insertion","stop_gained","stop_retained_variant",
                         "start_lost","stop_lost","coding_sequence_variant","incomplete_terminal_codon_variant")
CS.df.coding.PP01 <- CS.df.vep %>% filter(most_severe_consequence %in% coding_consequences,PP>0.01) %>% arrange(desc(PP))

# Write tables
write.table(CS.df.coding.PP01,file="../output/VEP/coding_variants_PP01.tsv",
            quote = FALSE, sep = "\t", col.names = T, row.names = F)
write.table(CS.df.vep,file="../output/VEP/CS_PP001_VEP_annotated.tsv",
            quote = FALSE, sep = "\t", col.names = T, row.names = F)


# CADD vs. PP
# Bin by PP
bins = c(0.001,0.01,0.05,0.25,0.75,1)
merged <- CS.df.vep
merged$PP <- as.numeric(as.character(merged$PP))
merged$PPbin <- cut(merged$PP, bins,include.lowest=T)
merged$CADD_PHRED <- as.numeric(as.character(merged$CADD_PHRED))

ggplot(merged %>% filter(PP>0.01),aes(x=PP,y=CADD_PHRED)) +
  geom_point() + 
  scale_color_manual(values =  jdb_palette("brewer_spectra")[-5]) +
  pretty_plot() + L_border() + 
  labs(x="PP") 














