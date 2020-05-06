library(tidyverse)
library(BuenColors)

# Choose Color Palette
THE_PALETTE <- jdb_palette("solar_rojos")

# Import data
gchromvar <- fread("../output/gchromVAR/gchromVAR_abf_MPN_arraycovar_meta_finngen_r4_zscores.txt", stringsAsFactors = FALSE)

# Set up coordinates
cellCoordsDF <- data.frame(
  CellLabel = c("HSC", "MPP", "LMPP", "CLP", "GMP-A", "GMP-B", "GMP-C", "CMP", "MEP", "NK", "CD4", "CD8", "B", "pDC", "Mono", "mDC", "Ery", "Mega"),
          x = c( 0,     0,      -5,    -5,      0,        -2,    2,       5,     7,    -10,   -8,    -6,   -4,  -2,     2,     4,      8,     10), 
          y = c(10,     8,      7,     5,       6,        5,     5,       7,     5,     2,     2,     2,    2,   2,     2,     2,      2,     2)
)

#--------------------------
# gchromVAR plots
#--------------------------

makeCVplot <- function(plottrait){
  df <- gchromvar
  plotdf <- merge(cellCoordsDF, df[df$Trait == plottrait, ],
                  by.x = "CellLabel", by.y = "Celltype")
  plotdf$zscore <- as.numeric(plotdf$zscore)
  plotdf$pvalue <-  pnorm(plotdf$zscore, lower.tail = FALSE)

  p1 <- ggplot(plotdf, aes(x = x, y = y, color =  -log10(pvalue))) + 
    geom_point(size = 11) + pretty_plot() +
    geom_text(aes(label=CellLabel),hjust=0.5, vjust=3) + 
    scale_color_gradientn(colors = THE_PALETTE, name = "-log10(pvalue)") +
    scale_y_continuous(limits = c(0, 11)) + ggtitle(plottrait)
  
  ggsave(p1, filename = paste0("../output/gchromVAR/rawPDFs/", plottrait, ".pdf"),
         height = 8, width = 10)
  return(plottrait)
}

plot_out <-  makeCVplot("MPN_arraycovar_meta_finngen_r4")

