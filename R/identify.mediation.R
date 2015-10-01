#' Identify Dot in Mediation Plot
#'
#' Enables user to click on a plot and get a correcsponting
#' line of the annotation data frame.
#'
#' @param medresults mediation object
#' @param label.col column of medresults to be plotted

#' @export
#' @examples
#' data(Tmem68)
#' med <- mediation.scan(target=Tmem68$target,
#'                       mediator=Tmem68$mediator,
#'                       annotation=Tmem68$annotation,
#'                       covar=Tmem68$covar,
#'                       qtl.geno=Tmem68$qtl.geno,
#'                       method="double-lod-diff")
#' plot(med, main="double-lod-diff")
#' identify(med)

identify.mediation <- function(medresults, label.col="Associated.Gene.Name"){
  names(medresults) = toupper(names(medresults))

  y = medresults$LOD
  x = gmb.coordinates(medresults$CHR, medresults$POS)

  label.col <- toupper(label.col)
  if (label.col %in% names(medresults)) labels <- medresults[,label.col] else labels <- medresults[,1]

  row <- identify(x=x, y=y, label=labels, n=1)
  med[row,]
}
