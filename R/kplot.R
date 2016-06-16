#' Karl Broman's interactive plot of conditioned LOD scores
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
#' kplot(med, symbol.col="symbol")
#'
#' # save the plot as HTML page
#' pl <- kplot(med, symbol.col="symbol")
#' tmpdir <- tempdir()
#' htmlwidgets::saveWidget(pl, paste0(tmpdir, "/Tmem68.html"),
#'                         selfcontained=FALSE, libdir = "_assets")
#'                         
#' @export

kplot <- function(med, symbol.col="symbol", chrlen=mouse.chrlen, ...){

  names(med) = toupper(names(med))
  symbol.col =  toupper(symbol.col)

  #Check input
  stopifnot(c("CHR", "POS","LOD") %in% names(med))
  if (!require(qtlcharts, quietly=TRUE)) stop("Requires 'qtlcharts' package, use devtools::install_github('kbroman/qtlcharts') to get it.")

  if (!("GMB" %in% names(med)))
    med$GMB <- gmb.coordinates(med$CHR, med$POS, chrlen=chrlen)
  if (symbol.col %in% names(med)) symbols <- med[,symbol.col] else symbols <- med[,1]

  # coloring of chromosomes
  max.chr <- max(as.numeric(med$CHR[grep("[0-9]+", med$CHR)])) # number of autosomes
  unique.chr <- levels(factor(factor(med$CHR, levels=c(1:max.chr, "X", "Y", "M")))) # all chromosomes, ordered 1..max.chr,X,Y,M
  chrcolor <- rep(6, length(unique.chr))
  chrcolor[1:length(unique.chr) %% 2 == 0] <- 3
  chrcolor <- factor(chrcolor, levels=1:6)
  names(chrcolor) <- unique.chr

  iplot(x=med$GMB, y=med$LOD,
        indID=symbols, group = chrcolor[med$CHR], 
        chartOpts = list(xlab="Genomic position (Mb)", ylab="LOD"),  ...)
}
