#' Karl Broman's Interactive Plot of conditioned LOD scores
#'
#'
#' Plot LOD statistics calculated by \code{mediation.scan} against
#' genomic positions. Display a gene name when hovering a mouse over.
#'
#'
#' @param medresults mediation object
#' @param symbol.col column with gene names
#' @param chrlen length of chromosomes as in \code{org.Mm.egCHRLENGTHS}

#' @export
#' @seealso \code{\link{plot.mediation}}, \code{\link{mediation.scan}}
#' @examples
#' data(Tmem68)
#' med <- mediation.scan(target=Tmem68$target,
#'                       mediator=Tmem68$mediator,
#'                       annotation=Tmem68$annotation,
#'                       covar=Tmem68$covar,
#'                       qtl.geno=Tmem68$qtl.geno)
#' kplot(med, symbol.col="Associated.Gene.Name")
#'
#' # save the plot as HTML page
#' pl <- kplot(med, symbol.col="Associated.Gene.Name")
#' tmpdir <- tempdir()
#' htmlwidgets::saveWidget(pl, paste0(tmpdir, "/Tmem68.html"),
#'                         selfcontained=FALSE, libdir = "_assets")

kplot <- function(medresults, symbol.col="Associated.Gene.Name", chrlen=mouse.chrlen, ...){

  names(medresults) = toupper(names(medresults))
  symbol.col =  toupper(symbol.col)

  #Check input
  stopifnot(c("CHR", "POS","LOD") %in% names(medresults))
  if (!require(qtlcharts, quietly=TRUE)) stop("Requires 'qtlcharts' package, use devtools::install_github('kbroman/qtlcharts') to get it.")

  if (!("GMB" %in% names(medresults)))
    medresults$GMB <- gmb.coordinates(medresults$CHR, medresults$POS, chrlen=chrlen)
  if (symbol.col %in% names(medresults)) symbols <- medresults[,symbol.col] else symbols <- medresults[,1]

  # coloring of chromosomes
  max.chr <- max(as.numeric(medresults$CHR[grep("[0-9]+", medresults$CHR)])) # number of autosomes
  unique.chr <- levels(factor(factor(medresults$CHR, levels=c(1:max.chr, "X", "Y", "M")))) # all chromosomes, ordered 1..max.chr,X,Y,M
  chrcolor <- rep(6, length(unique.chr))
  chrcolor[1:length(unique.chr) %% 2 == 0] <- 3
  chrcolor <- factor(chrcolor, levels=1:6)
  names(chrcolor) <- unique.chr

  iplot(x=medresults$GMB, y=medresults$LOD,
        indID=symbols, group = chrcolor[medresults$CHR], ...)
}
