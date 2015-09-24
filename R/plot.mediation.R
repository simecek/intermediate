#' Plots conditioned LOD scores against genome position
#'
#' For a given QTL haplotype probabilities `Q` and target `T`, the function tries to add all
#' possible mediators `M` and calculates LOD statistic. The low LOD value indicates `Q` and
#' `T` are conditionally independent given `M`, i.e. M is a mediator of causal relationship
#' for `Q` to T.
#'
#' @param medresults
#' @param chrlen length of chromosomes as in \code{org.Mm.egCHRLENGTHS}
#' @author Steven C. Munger

#' @export

plot.mediation = function(medresults, chrlen = mouse.chrlen, ...){

  #Check input
  stopifnot(c("CHR", "POS","LOD") %in% toupper(names(medresults)))

  names(medresults) = toupper(names(medresults))

  #Find unique chromosomes in results input, then pull out their lengths and calculate midpoints for plotting
  medresults$CHR = sub("X","20",medresults$CHR)
  medresults$CHR = sub("Y","21",medresults$CHR)
  medresults$CHR = sub("M","22",medresults$CHR)
  toMatch <- unique(medresults$CHR)
  toMatch = as.numeric(toMatch)
  toMatch = toMatch[order(toMatch)]
  toMatch = as.character(toMatch)
  toMatch = sub("20","X",toMatch)
  toMatch = sub("21","Y",toMatch)
  toMatch = sub("22","M",toMatch)
  chrlen = chrlen[match(toMatch,names(chrlen))]
  medresults$CHR = sub("20","X",medresults$CHR)
  medresults$CHR = sub("21","Y",medresults$CHR)
  medresults$CHR = sub("22","M",medresults$CHR)
  chrlen.copy = chrlen  #Create copy of chrlen


  ### Create X-axis gmb values for mediators (e.g. proteins/RNA) genome positions
  mb  = medresults$POS
  if(max(mb) > 3000){
    mb = medresults$POS*1e-6
  }
  gmb = mb
  medresults$CHR = as.character(medresults$CHR)
  unique.chr = unique(medresults$CHR)
  unique.chr = sub("X","20",unique.chr)
  unique.chr = as.numeric(unique.chr)
  unique.chr = unique.chr[order(unique.chr)]
  unique.chr = sub("20","X",unique.chr)
  #chrlen = chrlen[names(chrlen) %in% unique.chr]
  if(max(chrlen) > 200) {
    chrlen = chrlen * 1e-6
  } # if(max(chrlen) > 200)
  chrlen = cumsum(chrlen)
  chrlen = c(0, chrlen)
  chrmid = chrlen[-length(chrlen)] + (diff(chrlen) * 0.5)
  # Add the preceeding chromosome lengths to each SNP position.
  for(c in 2:length(unique.chr)) {
    rows = which(medresults$CHR == unique.chr[c])
    gmb[rows] = gmb[rows] + chrlen[c]
  } # for(c)
  gmb.mediation = gmb
  xmin.mediation = min(gmb.mediation)-30
  xmax.mediation = max(gmb.mediation)+30
  gmb.mediation.extended = c(xmin.mediation,gmb.mediation,xmax.mediation) #Used to extend the x-axis of the plot a bit


  chrlen = chrlen.copy
  if(max(chrlen) > 200) {
    chrlen = chrlen * 1e-6
  } # if(max(chrlen) > 200)
  chrlen = cumsum(chrlen)
  chrlen = c(0, chrlen)
  chrmid = chrlen[-length(chrlen)] + (diff(chrlen) * 0.5)


  #### Create the plot
  par(font = 2, font.lab = 2, font.axis = 2, xaxs="i",las = 1, mar=c(3, 4, 3, 1) + 0.1)
  plot(gmb.mediation.extended, c(0,medresults$LOD,0), col = 0, ylim = c(0, max(medresults$LOD)*1.05),ylab = "Conditioned LOD", xaxt = "n", xlab = "")
  usr = par("usr")
  rect(chrlen[2 * 1:(length(chrlen) * 0.5) - 1], usr[3],
       chrlen[2 * 1:(length(chrlen) * 0.5)], usr[4],
       border = "grey60", lty=3, col = rgb(0.96,0.96,0.96))
  rect(usr[1], usr[3], usr[2], usr[4], border = "black")
  points(gmb.mediation, medresults$LOD, type = "p", pch=21,lty = 1,lwd=1, cex=0.9, col="firebrick4") #dark red points
  text(chrmid, 0.97 * usr[4], names(chrlen)[-1], cex=1)


  return("Finished!")
}
