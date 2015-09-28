#' Identify Dot in Mediation Plot
#'
#'
#' @param medresults mediation object
#' @param label.col column of medresults to be plotted

#' @export

identify.mediation <- function(medresults, label.col="Associated.Gene.Name"){
  names(medresults) = toupper(names(medresults))

  y = medresults$LOD
  x = gmb.coordinates(medresults$CHR, medresults$POS)

  label.col <- toupper(label.col)
  if (label.col %in% names(medresults)) labels <- medresults[,label.col] else labels <- medresults[,1]

  row <- identify(x=x, y=y, label=labels, n=1)
  med[row,]
}
