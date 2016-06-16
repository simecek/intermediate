#' Identify dot in mediation plot
#'
#' Enable user to identify a dot on a plot produced by \code{\link{plot.mediation}}.
#'
#' @param med A mediation object
#' @param label.col A column name of \code{med} to be used to be plotted
#' 
#' @return A row from \code{med} corresponding to the point nearest to the click.
#' 
#' @seealso \code{\link{plot.mediation}}
#' 
#' @examples
#' data(Tmem68)
#' med <- mediation.scan(target = Tmem68$target,
#'                       mediator = Tmem68$mediator,
#'                       annotation = Tmem68$annotation,
#'                       covar = Tmem68$covar,
#'                       qtl.geno = Tmem68$qtl.geno)
#' plot(med)
#' identify(med)
#' @export

identify.mediation <- function(med, label.col="symbol"){
  names(med) = toupper(names(med))

  y = med$LOD
  x = gmb.coordinates(med$CHR, med$POS)

  label.col <- toupper(label.col)
  if (label.col %in% names(med)) {
    labels <- med[,label.col] 
  }  else {
    labels <- med[,1]
  }  

  row <- identify(x=x, y=y, label=labels, n=1)
  med[row,]
}
