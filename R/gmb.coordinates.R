#' Convert genomic coordinates
#'
#' Convert genomic coordinates from (chr, mb) representation to 
#' cummulative (gmb) representation. 
#' 
#' A cummulative (gmb) form could be useful for example as 
#' x-coordinate of a genomic plot.
#'
#' @param chr A character vector of chromosomes
#' @param pos A numeric vector of positions (b or Mb)
#' @param chrlen A vector with chromosome lengths (as in \code{org.Mm.egCHRLENGTHS})
#' 
#' @return A numeric vector with cummulative (gmb) genomic coordinates.
#' 
#' @seealso \code{\link{plot.mediation}}
#' @examples
#' gmb <- gmb.coordinates(Tmem68$annotation$chr, Tmem68$annotation$pos)
#' @export

gmb.coordinates <- function(chr, pos, chrlen = mouse.chrlen) {

  # chr must be numeric or X or Y or M
  chr <- as.character(chr)
  stopifnot(grepl("[0-9]+", chr) | chr %in% c("X", "Y", "M") )

  # length of all chromosomes must be given
  stopifnot(chr %in% names(chrlen))

  # pos must be finite, non-negative, numberic
  stopifnot(is.finite(pos) & is.numeric(pos) & pos>=0)

  # length of 'chr' vector equals length of 'pos' vector
  stopifnot(length(chr) == length(pos))

  # number of autosomes
  max.chr <- max(as.numeric(chr[grep("[0-9]+", chr)]))

  # all chromosomes, ordered 1..max.chr,X,Y,M
  unique.chr <- levels(factor(factor(chr, levels=c(1:max.chr, "X", "Y", "M"))))

  # chrlen ordered as unique.chr
  chrlen = chrlen[unique.chr]

  # if pos > MAXPOS, then pos is NOT in Mb and needs conversion
  MAXPOS = 3000
  if (max(chrlen) > MAXPOS) chrlen <- chrlen / 10^6
  if (max(pos) > MAXPOS) pos <- pos / 10^6

  # cumulative chr. lengths (=shift)
  chrcumsum <- c(0, cumsum(chrlen))
  names(chrcumsum) <- c(names(chrlen), "End") # shift by 1

  gmb <- pos + chrcumsum[chr]
  return(gmb)
}
