#' Interactive plot of conditioned LOD scores
#'
#' Plot LOD statistics calculated by \code{mediation.scan} against
#' genomic positions. Display a gene name when hovering a mouse over.
#'
#' @return Htmlwidget \code{iplot}
#'
#' @param med A mediation object
#' @param symbol.col A column of \code{med} with gene names
#' @param chrlen A vector with chromosome lengths (as in \code{org.Mm.egCHRLENGTHS})

#' @seealso \code{\link{plot.mediation}}, \code{\link{mediation.scan}}
#' 
#' @examples
#' data(Tmem68)
#' med <- mediation.scan(target = Tmem68$target,
#'                       mediator = Tmem68$mediator,
#'                       annotation = Tmem68$annotation,
#'                       covar = Tmem68$covar,
#'                       qtl.geno = Tmem68$qtl.geno)
#' 
#' theme_set(theme_bw())
#' kplot(med, symbol.col="symbol")
#'
#' # save the plot as HTML page
#' gp <- kplot(med, symbol.col="symbol") 
#' pl <- ggplotly(gp, tooltip="text")
#' pl
#' tmpdir <- tempdir()
#' htmlwidgets::saveWidget(pl, paste0(tmpdir, "/Tmem68.html"),
#'                         selfcontained=TRUE)
#'                         
#' @export

kplot <- function(med, symbol.col="symbol", chrlen=mouse.chrlen, ...){
  
  names(med) = toupper(names(med))
  symbol.col =  toupper(symbol.col)
  
  #Check input
  stopifnot(c("CHR", "POS","LOD") %in% names(med))
  
  if (!("GMB" %in% names(med)))
    med$GMB <- gmb.coordinates(med$CHR, med$POS, chrlen=chrlen)
  if (symbol.col %in% names(med)) symbols <- med[,symbol.col] else symbols <- med[,1]
  
  # coloring of chromosomes
  max.chr <- max(as.numeric(med$CHR[grep("[0-9]+", med$CHR)])) # number of autosomes
  unique.chr <- levels(factor(factor(med$CHR, levels=c(1:max.chr, "X", "Y", "M")))) # all chromosomes, ordered 1..max.chr,X,Y,M
  chrcolor <- rep(1, length(unique.chr))
  chrcolor[1:length(unique.chr) %% 2 == 0] <- 2
  chrcolor <- factor(chrcolor, levels=1:6)
  names(chrcolor) <- unique.chr
  
  med$color = chrcolor[med$CHR]
  
  # calculates breakpoints between chromosomes
  get_chr_breaks <- function(chr, gmb) {
    tmp <- data.frame(chr, gmb)
    tmp2 <- tmp %>% group_by(chr) %>% summarize(maxgmb = max(gmb))
    gmb.breakpoints <- c(0, sort(tmp2$maxgmb))
    return(sort(gmb.breakpoints))
  }
  
  # calculates where to plot chromosome name
  get_chr_middle_points <- function(chr, gmb) {
    tmp <- data.frame(chr, gmb)
    tmp2 <- tmp %>% group_by(chr) %>% summarize(maxgmb = min(gmb)/2 + max(gmb)/2)
    gmb.breakpoints <- sort(tmp2$maxgmb)
    return(sort(gmb.breakpoints))
  }
  
  chrs <- c(as.character(1:19), "X", "Y", "M")
  chrs <- chrs[chrs %in% unique(med$CHR)]
  gene.breaks <- get_chr_breaks(med$CHR, med$GMB)
  gene.mid <- get_chr_middle_points(med$CHR, med$GMB)
  
  # to be corrected !!!
  line.style <- theme_bw()$panel.grid.major
  
  ggplot(med, aes(x=GMB, y=LOD, color=color, text=SYMBOL)) +
    geom_point() +
    scale_x_continuous(breaks=gene.mid, labels = chrs, minor_breaks=gene.breaks, expand=c(0,0)) +
    theme(legend.position="none") +
    theme(panel.grid.major.x = element_blank()) + 
    theme(panel.grid.minor.x = line.style) +
    xlab('Position')
  
}