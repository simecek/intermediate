#' Plots conditioned LOD scores against genomic positions
#'
#' Plot LOD statistics calculated by \code{mediation.scan} against genomic positions.
#'
#' @param med A mediation object
#' @param col A color of points
#' @param chrlen A vector with chromosome lengths (as in \code{org.Mm.egCHRLENGTHS})
#' 
#' @seealso \code{\link{identify.mediation}}, \code{\link{kplot}}
#' 
#' @examples
#' data(Tmem68)
#' med <- mediation.scan(target=Tmem68$target,
#'                       mediator=Tmem68$mediator,
#'                       annotation=Tmem68$annotation,
#'                       covar=Tmem68$covar,
#'                       qtl.geno=Tmem68$qtl.geno)
#' plot(med)
#' @export

plot.mediation <- function(med, 
                           col="firebrick4", 
                           chrlen=mouse.chrlen, ...){

  names(med) = toupper(names(med))

  #Check input
  stopifnot(c("CHR", "POS","LOD") %in% names(med))

  if (!("GMB" %in% names(med)))
    med$GMB <- gmb.coordinates(med$CHR, med$POS, chrlen=chrlen)

  # reorganize chr-lengths as in gmb.coordinates
  max.chr <- max(as.numeric(med$CHR[grep("[0-9]+", med$CHR)]))
  unique.chr <- levels(factor(factor(med$CHR, levels=c(1:max.chr, "X", "Y", "M"))))
  chrlen = chrlen[unique.chr]

  # if pos > MAXPOS, then pos is NOT in Mb and needs conversion
  MAXPOS = 3000
  if (max(chrlen) > MAXPOS) chrlen <- chrlen / 10^6

  ### Create X-axis gmb values for mediators (e.g. proteins/RNA) genome positions
  gmb.mediation.extended = c(min(med$GMB)-30,med$GMB,max(med$GMB)+30) #Used to extend the x-axis of the plot a bit
  chrlen = c(0, cumsum(chrlen))
  chrmid = chrlen[-length(chrlen)] + (diff(chrlen) * 0.5)

  #### Create the plot
  par(font = 2, font.lab = 2, font.axis = 2, xaxs="i",las = 1, mar=c(3, 4, 3, 1) + 0.1)
  plot(gmb.mediation.extended, c(0,med$LOD,0), col = 0, ylim = c(0, max(med$LOD)*1.05),ylab = "Conditioned LOD", xaxt = "n", xlab = "", ...)
  usr = par("usr")
  rect(chrlen[2 * 1:(length(chrlen) * 0.5)-1], usr[3],
       chrlen[2 * 1:(length(chrlen) * 0.5)], usr[4],
       border = "grey60", lty=3, col = rgb(0.96,0.96,0.96))
  rect(usr[1], usr[3], usr[2], usr[4], border = "black")
  points(med$GMB, med$LOD, type = "p", pch=21,lty = 1,lwd=1, cex=0.9, col=col) #dark red points
  text(chrmid, 0.97 * usr[4], names(chrlen)[-1], cex=1)
}
