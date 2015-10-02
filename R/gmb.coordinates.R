#' For chromosome and position, calculates a genomic coordinate
#'
#' Convert (chr, mb) coordinates to (gmb) coordinates
#'
#'
#' @param chr chromosome
#' @param pos position (b or Mb)
#' @param chrlen lengths of chromosomes as in \code{org.Mm.egCHRLENGTHS}

#' @export
#' @seealso \code{\link{plot.mediation}}
#' @examples
#' gmb <- gmb.coordinates(Tmem68$annotation$Chr, Tmem68$annotation$Pos)

gmb.coordinates <- function(chr, pos, chrlen=mouse.chrlen) {

  # chr must be numeric or X or Y or M
  chr <- as.character(chr)
  stopifnot(grepl("[0-9]+", chr) | chr %in% c("X", "Y", "M") )

  # length of all chrs must be given
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

  # middle chr. position
  chrmid <- chrcumsum[-length(chrcumsum)] + (diff(chrcumsum) * 0.5)

  gmb <- pos + chrcumsum[chr]
  return(gmb)
}
