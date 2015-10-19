#' Mediation Scan
#'
#' For a given QTL haplotype probabilities \code{qtl.geno} and target \code{target},
#' the function sequentially tries to add each column \code{m} of \code{mediator} matrix as a covariate
#' and calculates LOD statistic. The low LOD value indicates \code{qtl.geno} and
#' \code{target} are conditionally independent given \code{m},
#' i.e. \code{m} is a mediator of causal relationship from \code{qtl.geno} to \code{target}.
#'
#'
#' @param target a numeric vector with gene/protein expression
#' @param mediator a matrix, each column is one gene/protein's expression
#' @param annotation a data frame with mediators' annotation, must include columns "chr" and "pos"
#' @param qtl.geno a matrix, haplotype probabilities at QTL we try to mediate
#' @param covar a matrix with additive covariates
#' @param method a method to handle missing cases
#' @param verbose if TRUE display information about the progress

#' @export
#' @seealso \code{\link{plot.mediation}}, \code{\link{kplot}}

#' @examples
#' data(Tmem68)
#' med <- mediation.scan(target=Tmem68$target,
#'                       mediator=Tmem68$mediator,
#'                       annotation=Tmem68$annotation,
#'                       covar=Tmem68$covar,
#'                       qtl.geno=Tmem68$qtl.geno,
#'                       method="double-lod-diff")
#' plot(med, main="double-lod-diff")

mediation.scan <- function(target, mediator, annotation, qtl.geno, covar=NULL,
                           method=c("double-lod-diff", "ignore", "lod-diff", "lod-ratio"), verbose=TRUE) {

  # calculates log10-Likelihood of linear model y ~ 1 + X
  LL <- function(y, X) {
    -length(y)/2*log10(sum(qr.resid(qr(cbind(X,1)),y)^2))
  }

  # check input
  stopifnot(NROW(target) == NROW(mediator))
  stopifnot(NROW(annotation) == NCOL(mediator))
  stopifnot(NROW(qtl.geno) == NROW(target))
  stopifnot(is.null(covar) | NROW(target) == NROW(covar))
  stopifnot(!any(is.na(covar)))
  stopifnot(!any(is.na(qtl.geno)))
  stopifnot(all(is.numeric(target)))
  stopifnot(all(is.numeric(mediator)))
  stopifnot(all(is.numeric(qtl.geno)))
  stopifnot(all(is.numeric(covar)))
  stopifnot(c("CHR", "POS") %in% toupper(names(annotation)))
  method = match.arg(method)

  # data preparation
  mediator <- cbind(mediator) # to ensure 'mediator' is a matrix
  N <- ncol(mediator) # number of points to scan
  if (is.null(covar)) covar <- cbind(rep(1, N)) # if no covariates, use just intercept
  LOD <- rep(NA, N) # prepare output

  if (method == "double-lod-diff") {
    no.na <- !is.na(target)
    LOD0 <- LL(target[no.na], cbind(covar, qtl.geno)[no.na,]) - LL(target[no.na], covar[no.na,])
  }

  # for-loop comparing M0: target~covar+mediator[,i] vs M1: target~covar+mediator[,i]+qtl.geno
  for (i in 1:N) {
    if (verbose & i %% 1000 == 0) print(i)
    no.na <- !is.na(target) & !is.na(mediator[,i])
    loglik0 <- LL(target[no.na], cbind(covar[no.na,], mediator[no.na,i]))
    loglik1 <- LL(target[no.na], cbind(covar[no.na,], mediator[no.na,i], qtl.geno[no.na,]))

    if (method == "ignore" | (method == "double-lod-diff" & all(no.na))) {
      # "double-lod-diff" for no missing observation is identical to "ignore"
      LOD[i] <- loglik1 - loglik0
    } else {
      loglik2 <- LL(target[no.na], covar[no.na,])
      loglik3 <- LL(target[no.na], cbind(covar[no.na,], qtl.geno[no.na,]))

      if (method == "lod-diff") {
        LOD[i] <- loglik3 - loglik2 - (loglik1-loglik0)
      } else if (method == "double-lod-diff") {
        LOD[i] <- LOD0 - (loglik3 - loglik2 - (loglik1-loglik0))
      } else if (method == "lod-ratio") {
        LOD[i] <- (10^loglik1-10^loglik0) / (10^loglik3 - 10^loglik2)
      }
    }
  }

  output <- annotation
  output$LOD <- LOD
  class(output) <- c("mediation", "data.frame")
  return(output)
}
