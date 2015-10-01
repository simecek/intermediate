#' Plots conditioned LOD scores against genome position
#'
#' For a given QTL haplotype probabilities `Q` and target `T`, the function tries to add all
#' possible mediators `M` and calculates LOD statistic. The low LOD value indicates `Q` and
#' `T` are conditionally independent given `M`, i.e. M is a mediator of causal relationship
#' for `Q` to T.
#'
#' @param medresults mediation object
#' @param verbose if TRUE then prints "Finished!" when done
#' @param chrlen length of chromosomes as in \code{org.Mm.egCHRLENGTHS}

#' @author Steven C. Munger

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

plot.mediation <- function(medresults, col="firebrick4", chrlen=mouse.chrlen, verbose=FALSE, ...){

  names(medresults) = toupper(names(medresults))

  #Check input
  stopifnot(c("CHR", "POS","LOD") %in% names(medresults))

  if (!("GMB" %in% names(medresults)))
    medresults$GMB <- gmb.coordinates(medresults$CHR, medresults$POS, chrlen=chrlen)

  # reorganize chr-lengths as in gmb.coordinates
  max.chr <- max(as.numeric(medresults$CHR[grep("[0-9]+", medresults$CHR)]))
  unique.chr <- levels(factor(factor(medresults$CHR, levels=c(1:max.chr, "X", "Y", "M"))))
  chrlen = chrlen[unique.chr]

  # if pos > MAXPOS, then pos is NOT in Mb and needs conversion
  MAXPOS = 3000
  if (max(chrlen) > MAXPOS) chrlen <- chrlen / 10^6

  ### Create X-axis gmb values for mediators (e.g. proteins/RNA) genome positions
  gmb.mediation.extended = c(min(medresults$GMB)-30,medresults$GMB,max(medresults$GMB)+30) #Used to extend the x-axis of the plot a bit
  chrlen = c(0, cumsum(chrlen))
  chrmid = chrlen[-length(chrlen)] + (diff(chrlen) * 0.5)

  #### Create the plot
  par(font = 2, font.lab = 2, font.axis = 2, xaxs="i",las = 1, mar=c(3, 4, 3, 1) + 0.1)
  plot(gmb.mediation.extended, c(0,medresults$LOD,0), col = 0, ylim = c(0, max(medresults$LOD)*1.05),ylab = "Conditioned LOD", xaxt = "n", xlab = "", ...)
  usr = par("usr")
  rect(chrlen[2 * 1:(length(chrlen) * 0.5)-1], usr[3],
       chrlen[2 * 1:(length(chrlen) * 0.5)], usr[4],
       border = "grey60", lty=3, col = rgb(0.96,0.96,0.96))
  rect(usr[1], usr[3], usr[2], usr[4], border = "black")
  points(medresults$GMB, medresults$LOD, type = "p", pch=21,lty = 1,lwd=1, cex=0.9, col=col) #dark red points
  text(chrmid, 0.97 * usr[4], names(chrlen)[-1], cex=1)

  if (verbose) return("Finished!")
}
