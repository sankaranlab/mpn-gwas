library(gtools)
library(ggrastr)

lambda.gc <- function(p,plot=FALSE,col = palette()[4], lcol = palette()[2], ...) {
  ## for any vector p containing p-values
  ## Convert p-value to Chi-sq and calculate genomic controls
  ## Modify from gcontrol2() from library("gap") version 1.1.3
  ## Bhoom Suktitipat
  ## 08/01/2012
  p <- p[!is.na(p)]
  n <- length(p)
  x2obs <- qchisq(p, 1, lower.tail = FALSE)
  x2exp <- qchisq(1:n/n, 1, lower.tail = FALSE)
  lambda <- median(x2obs)/median(x2exp)
  if (plot) qqunif(p, col = col, lcol = lcol, ...)
  return(lambda)
}

#' Creates a Q-Q plot
#' 
#' Creates a quantile-quantile plot from p-values from a GWAS study.
#' 
#' @param pvector A numeric vector of p-values.
#' @param ... Other arguments passed to \code{plot()}
#' 
#' @return A Q-Q plot.
#' 
#' @keywords visualization qq qqplot
#' 
#' @importFrom stats ppoints
#' @import utils
#' @import graphics
#' 
#' @examples
#' qq(gwasResults$P)
#' 
#' @export

qq = function(pvector, ...) {
  
  # Check for sensible input
  if (!is.numeric(pvector)) stop("Input must be numeric.")
  
  # limit to not missing, not nan, not null, not infinite, between 0 and 1
  pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector<1 & pvector>0]
  
  # Observed and expected
  o = -log10(sort(pvector,decreasing=FALSE))
  e = -log10( ppoints(length(pvector) ))
  
  
  #     # The old way
  #     plot(e, o, pch=20, 
  #          xlab=expression(Expected~~-log[10](italic(p))), 
  #          ylab=expression(Observed~~-log[10](italic(p))), 
  #          ...)
  
  # The new way to initialize the plot.
  ## See http://stackoverflow.com/q/23922130/654296
  ## First, define your default arguments
  def_args <- list(pch=20, xlim=c(0, max(e)), ylim=c(0, max(o)), 
                   xlab=expression(Expected~~-log[10](italic(p))), 
                   ylab=expression(Observed~~-log[10](italic(p)))
  )
  ## Next, get a list of ... arguments
  #dotargs <- as.list(match.call())[-1L]
  dotargs <- list(...)
  ## And call the plot function passing NA, your ... arguments, and the default
  ## arguments that were not defined in the ... arguments.
  tryCatch(do.call("plot", c(list(x=e, y=o), def_args[!names(def_args) %in% names(dotargs)], dotargs)), warn=stop)
  
  # Add diagonal
  abline(0,1,col="red")
  
}


#' Creates a Manhattan Plot.
#'
#' Creates a Manhattan Plot as a ggplot layer.
#'
#' The function creates a manhattan plot as a ggplot layer. The output can be stored in a variable,
#' and additional layers can be added. See \code{\link{ggmanHighlight}}, \code{\link{ggmanZoom}} and
#'  \code{\link{ggmanLabel}}.
#' 
#' @param gwas A data frame with the gwas results
#' @param clumps Optional argument; takes an object of class 'ggclumps' containing
#' the SNP clumps, see \code{\link{ggmanClumps}}
#' @param clumps.colour colour of the clumps
#' @param snp Name of the column containing SNP identifiers; default is 'snp'
#' @param bp Name of the column containing the basepair positions; default is 'bp'
#' @param chrom Name of the column containing the chromosome identifers; default is 'chrom'
#' @param pvalue Name of the column containing the p values; default is 'pvalue'
#' @param sigLine Threshold for the genome wide significant line in -log10 scale; default is 8; specify NA to remove the line.
#' @param lineColour colour of the genomewide line; default is 'red'
#' @param pointSize Size of the points in the plot; default is 0.1
#' @param ymin Starting point of y axis; default is 0
#' @param ymax Ending point of y axis
#' @param logTransform if TRUE, P value is -log10 transformed; default is TRUE; Specify FALSE
#' when plotting values other p values, such as zscore or beta
#' @param invert if TRUE, an inverted Manhattan Plot will be created. The p values of the variants with (or < 1 or beta < 0) will positive log10-transformed, which will result in negative values. 
#' @param invert.method whether inversion should be based on odds ratio or beta. possible values: 'or' or 'beta'
#' @param invert.var name of the column in the gwas data.frame containing the or or beta
#' @param relative.positions if TRUE, the x axis points will be calculated in proportion to the basepair positions. So, the gaps in the genome with no genotypes will be reflected in the plot(requires more computation, hence more time to plot). If FALSE, the SNPs are ascendingly sorted chromosome-wise. The default value is FALSE. 
#' @param xlabel X-axis label
#' @param ylabel Y-axis label
#' @param title plot title
#' @param legend.title title of legend; this argument applies only when the clumps are plotted and highlighted in groups.
#' @param clumps.label.type type of label; either 'text' or 'label'; default is 'label'
#' @param legend.remove if TRUE, the  legend will be removed; default is FALSE
#' @param ... other arguments to pass to \code{\link[ggplot2]{geom_point}};
#' Note: do not use \emph{size} to specify size of the points, instead use \emph{pointSize}
#'
#' @import ggplot2
#' @import ggrepel
#' 
#' @importFrom gtools mixedorder
#' 
#'
#' @return A Manhattan Plot
#'
#' 
#' @examples
#'
#'#simple Manhattan Plot
#' ggman(toy.gwas, snp = "snp", bp = "bp", chrom = "chrom",
#' pvalue = "pvalue")
#'
#' #enable relative positioning
#' ggman(toy.gwas, snp = "snp", bp = "bp", chrom = "chrom",
#' pvalue = "pvalue",relative.positions = TRUE)
#'
#' #plot odds ratio
#' ggman(toy.gwas, snp = "snp", bp = "bp", chrom = "chrom",
#' pvalue = "or",logTransform = FALSE, ymax = 3)
#' 
#' #plot beta
#' ggman(toy.gwas, snp = "snp", bp = "bp", chrom = "chrom", pvalue = "beta",
#' logTransform = FALSE, ymin = -2, ymax = 2)
#'
#' #inverted Manhattan plot 
#' ggman(toy.gwas, snp = "snp", bp = "bp", chrom = "chrom", pvalue = "pvalue",
#' invert = TRUE, invert.method = 'or', invert.var = "or")
#' 
#'
#' 
#' 
#'
#' 
#'
#'
#' @export
ggman <- function(gwas,
                  clumps = NA,
                  clumps.colour = "blue",
                  snp = NA, bp = NA, chrom = NA, pvalue = NA,index=NA,
                  sigLine = 8,
                  lineColour = "red",
                  pointSize = 0.2,
                  ymin = NA, ymax = NA,
                  logTransform = TRUE,
                  invert = FALSE, invert.method='or',invert.var='or',
                  relative.positions = FALSE,
                  xlabel = "chromosome", ylabel = "-log10(P value)", title = "Manhattan Plot",
                  legend.title = "legend", clumps.label.type = 'label', legend.remove = FALSE,
                  odd.xlabels=FALSE,...) {
  
  ##define global variables to escape R CMD check
  
  
  dfmnames <- names(gwas)
  ## chrom input
  if(is.na(chrom)){
    chrom <- search.names(c("chr","chrom","chromosome"), dfmnames)
    if(is.null(chrom)){
      stop("Couldn't find the chromosome column.
           Specify the name of the column with chromosome ids")
    }
    }
  
  if(is.na(snp)){
    snp <- search.names(c("snp","rsid"), dfmnames)
    if(is.null(snp)){
      stop("Couldn't find the snp column.
           Specify the name of the column with snp ids")
    }
    }
  
  if(is.na(bp)){
    bp <- search.names(c("bp","basepair","position","start"),dfmnames)
    if(is.null(bp)){
      stop("Couldn't find the bp column.
           Specify the name of the column with bp ids")
    }
    }
  if(is.na(pvalue)){
    pvalue <- search.names(c("p","pval","pvalue"), dfmnames)
    if(is.null(pvalue)){
      stop("Couldn't find the pvalue column.
           Specify the name of the column with pvalues")
    }
    }
  

  dfm <- as.data.frame(gwas)
  dfm$chrom <- dfm[,chrom]
  dfm$bp <- as.numeric(as.character(dfm[,bp]))
  dfm$pvalue <- as.numeric(as.character(dfm[,pvalue]))
  dfm$snp <- dfm[,snp]
  dfm$chrom <- as.character(dfm$chrom)
  
  if (is.na(index)){
    ##add index
    dfm <- dfm[order(dfm$bp),]
    dfm <- dfm[mixedorder(dfm$chrom),]
    dfm$index <- 1:nrow(dfm)
  } else {
    dfm$index <- dfm[,index]
    dfm <- dfm[order(dfm$index),]
    dfm <- dfm[mixedorder(dfm$chrom),]
  }
  
  ##find the number of chromosomes
  chrtable <- data.frame(table(dfm$chrom))
  chrtable$Var1 <- as.character(chrtable$Var1)
  chrtable <- chrtable[mixedorder(chrtable$Var1),]
  ##group odd and even chromosomes
  oddchrom <- as.character(chrtable$Var1[seq(1,nrow(chrtable),2)])
  dfm$chrom_alt <- replace(dfm$chrom, dfm$chrom %in% oddchrom, 0)
  dfm$chrom_alt <- replace(dfm$chrom_alt, dfm$chrom_alt != 0,1)
  
  if(logTransform){
    dfm$marker <- -log10(dfm$pvalue)
  } else {
    dfm$marker <- dfm$pvalue
  }
  
  if (is.na(ymax)){
    ymax <- max(-log10(dfm$pvalue)) + 1
  }
  
  if(relative.positions){
    relpos <- function(x,minbp,maxbp,nrows,startingpoint){
      actual.distance <- (x - minbp)
      relative.distance <- (actual.distance*nrows)/maxbp
      return(relative.distance + startingpoint)
    }
    dfm$chrom <- factor(dfm$chrom, levels = chrtable$Var1)
    dfm.split <- split(dfm, dfm$chrom)
    startingpoint = 1
    dfm.list <- lapply(dfm.split, function(x){         
      minbp <- as.numeric(min(x$bp))
      maxbp <- as.numeric(max(x$bp))
      nrows <- as.numeric(nrow(x))
      x$index <- relpos(x$bp,minbp,maxbp,nrows,startingpoint)
      startingpoint <<- max(x$index)+1
      return(x)
    })
    dfm <- do.call(rbind,dfm.list)
  }
  
  ##create x axis tick points
  dfmsplit <- split(dfm, dfm$chrom)
  xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    ##x$index[length(x$index)/2]
    return(x$index[midpoint])
  })
  
  if(odd.xlabels){
    ##sort the breaks base on names
    xbreaks <- xbreaks[mixedorder(names(xbreaks))]
    xbreaks <- xbreaks[seq(1,length(xbreaks),2)]
    print(xbreaks)
  }
  ##invert
  if(invert){
    if(invert.method == 'or'){
      dfm$or <- as.numeric(as.character(dfm[,invert.var]))
      dfm$sign <- with(dfm, replace(pvalue,or > 1, 1))
      dfm$sign <- with(dfm, replace(sign, sign != 1, -1))
    } else {
      dfm$beta <- as.numeric(as.character(dfm[,invert.var]))
      dfm$sign <- with(dfm, replace(beta, beta > 0, 1))
      dfm$sign <- with(dfm, replace(sign, sign != 1, -1))
    }
    dfm$marker <- with(dfm, -log10(pvalue) * sign)
    ##axis
    if(is.na(ymin)){
      ymin = ymax * -1
    }        
  }
  
  
  if (!is.na(clumps)[1]){
    if (!any(class(clumps) == "ggmanClumps")) {
      stop("clumps argument takes an object of class ggmanclumps;see ggmanClumps function")
    }
    environment(plot.clumps) <- environment()
    plot.clumps()        
  } else {
    p1 <- ggplot(dfm, aes(x = index,y = marker, colour = as.factor(chrom_alt))) +
      geom_point_rast(size=pointSize) +
      scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),
                         expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0), limits = c(ymin,ymax))+
      guides(colour = FALSE) +
      labs(x = xlabel, y = ylabel, title = title)
    
    if(!is.na(sigLine)){
      p1 <- p1 + geom_hline(aes(yintercept= as.numeric(sigLine)),
                            colour = lineColour, size = 0.25,linetype="dashed") +
        geom_hline(aes(yintercept= as.numeric(-log10(1e-6))),
                   colour = "blue", size = 0.25,linetype="dashed") 
    }
    if(invert){
      p1 <- p1 + geom_hline(aes(yintercept= as.numeric(sigLine) * -1),
                            colour = lineColour, size = 0.25) +
        geom_hline(aes(yintercept= 0),colour = "white", size = 0.25)
    }
    class(p1) <- append(class(p1), "ggman")
    return(p1)
  }
  
}

ggmanHighlight <- function(ggmanPlot, highlight,colour = "blue", ...){
  ##input checks
  environment(check.input.ggmanHighlight) <- environment()
  check.input.ggmanHighlight()
  
  dfm <- ggmanPlot[[1]]
  dfm <- dfm[dfm$snp %in% highlight,]
  if(nrow(dfm) == 0){
    stop("None of the markers in the input is present in the Manhattan plot layer")
  }
  ggmanPlot +
    scale_colour_grey(start = 0.5,end = 0.6) +
    geom_point_rast(data = dfm,colour= colour,...)
}