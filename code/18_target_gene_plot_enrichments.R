library(data.table)
library(tidyverse)
library(BuenColors)

# Read gene set enrichments
fuma <- fread("../output/target_genes/target_gene_FUMA_coding_variants_CHEK2_included.txt") %>% filter(Category == "GO_bp") %>% arrange(adjP)
fuma$GeneSet <- gsub("GO_","",fuma$GeneSet) %>% gsub("_", " ",.) %>% str_to_title()
fuma$GeneSet <- factor(fuma$GeneSet,levels=rev(fuma$GeneSet))

# Plot
p1 <- ggplot(head(fuma,10),aes(x=GeneSet,label=GeneSet,y=-log10(p)))+
  geom_bar(stat="identity",fill=jdb_palette("solar_extra")[4],position = position_stack(reverse = TRUE)) + 
  coord_flip() +
  theme_bw() +
  pretty_plot(fontsize = 8) +
  L_border()+
  scale_y_continuous(expand=c(0,0))+
  theme(legend.position="none",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  labs(x="",y="-log(pvalue)") +
  geom_text(color="black",
            position = position_stack(vjust = 0.02),hjust=0,size=4)
p1
cowplot::ggsave2(p1, file="../output/target_genes/fuma_GO_enrichments.pdf",width=3,height=2.5)

