#' @title Maximum Number of Mismatches
#'
#' @description Calculate the maximum expected number of mismatches for
#'   duplicate samples, parent-offspring pairs, and parent-parent-offspring
#'   trios.
#'
#' @details The thresholds for maximum number of mismatches calculated here aim
#'   to minimise false negatives, i.e. to minimise the chance that any true
#'   duplicates or true parent-offspring pairs are already excluded during the
#'   filtering steps where these \code{MaxMismatch} values are used.
#'   Consequently, there is a high probability of false positives, i.e. it is
#'   likely that some sample pairs with fewer mismatches than the
#'   \code{MaxMismatch} threshold, are in fact not duplicate samples or
#'   parent-offspring pairs. Use of these \code{MaxMismatch} thresholds is
#'   therefore only the first step of pedigree reconstruction by
#'   \code{\link{sequoia}}.
#'
#' @param Err estimated genotyping error rate, as a single number or 3x3 matrix
#'   (averaged value(s) across SNPs), or a vector with the same length as MAF,
#'   or a nSnp x 3 x 3 array. If a matrix, this should be the probability of
#'   observed genotype (columns) conditional on actual genotype (rows). Each row
#'  must therefore sum to 1. If an array, each 3x3 slice should abide this rule.
#' @param MAF vector with minor allele frequency at each SNP.
#' @param ErrFlavour function that takes \code{Err} as input, and returns a 3x3
#'   matrix of observed (columns) conditional on actual (rows) genotypes, or
#'   choose from inbuilt ones as used in sequoia 'version2.0', 'version1.3', or
#'   'version1.1'. Ignored if \code{Err} is a matrix. See \code{\link{ErrToM}}.
#' @param qntl quantile of binomial distribution to be used as the maximum, of
#'   individual-level probability. For a desired dataset-level probability
#'   quantile \eqn{Q}, use \code{qntl}\eqn{= Q^{(1/N)}}, where \eqn{N} is
#'   the number of individuals.
#' @param Return Either 'Counts' to return the threshold counts (default),
#'  or 'Probs' to return the mismatch probabilities from which these counts
#'  are calculated.
#'
#' @return A vector with three integers:
#'  \item{DUP}{Maximum number of differences between 2 samples from the
#'   same individual}
#'  \item{OH}{Maximum number of Opposing Homozygous SNPs between a true
#'    parent-offspring pair}
#'  \item{ME}{Maximum number of Mendelian Errors among a true parent-parent-
#'   offspring trio}.
#'
#' @seealso  \code{\link{SnpStats}}.
#'
#' @importFrom stats setNames
#'
#' @examples
#' CalcMaxMismatch(Err = 0.05, MAF = runif(n=100, min=0.3, max=0.5))
#' \dontrun{
#' CalcMaxMismatch(Err = 0.02, MAF = SnpStats(MyGenoMatrix, Plot=FALSE)[,"AF"])
#' }
#'
#' @export

CalcMaxMismatch <- function(Err, MAF, ErrFlavour = "version2.9", qntl=1-1e-5,
                            Return='Counts') {

  nSnp <- length(MAF)
  # ErrM: observed (columns) conditional on actual (rows)
  if (length(Err)==1 | (length(Err)==9 & all(dim(Err)==c(3,3)))) {
    ErrM <- ErrToM(Err, flavour = ErrFlavour, Return = "matrix")  # performs checks if Err is already matrix
    ErrA <- aperm(array(ErrM, dim=c(3,3, nSnp)), c(3,1,2))
  } else if (length(Err) == nSnp) {
    ErrA <- plyr::laply(Err, ErrToM, flavour = ErrFlavour, Return = "matrix")
  } else if (length(Err) == 9*nSnp & all(dim(Err) == c(nSnp,3,3))) {
    ErrA <- plyr::aaply(Err, ErrToM, flavour = ErrFlavour, Return = "matrix")  # check each slice
  } else {
    stop("'Err' must be a single value, or a vector with same length as 'MAF'")
  }

  if (!is.double(MAF) || any(MAF < 0 | MAF > 1)) {
    stop("Allele frequency vector 'MAF' must be between 0 and 1")
  }
  if (!is.double(qntl) || qntl < 0 | qntl > 1) {
    stop("'qntl' must be a number between 0 and 1")
  }


  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # arrays to calculate offspring-dam-sire mendelian error probability
  AKA2P <- array(dim=c(3,3,3))  # offspr - mother - father
  AKA2P[1,,] <- matrix(c(1,.5,0, .5,.25, 0,  0, 0,0), nrow=3)
  AKA2P[2,,] <- matrix(c(0,.5,1, .5, .5,.5,  1,.5,0), nrow=3)
  AKA2P[3,,] <- matrix(c(0, 0,0,  0,.25,.5,  0,.5,1), nrow=3)

  I.ME <- array(0, dim=c(3,3,3))
  I.ME[1,3, ] <- 1
  I.ME[1, ,3] <- 1
  I.ME[1,3,3] <- 2
  I.ME[3,1, ] <- 1
  I.ME[3, ,1] <- 1
  I.ME[3,1,1] <- 2
  I.ME[2,1,1] <- 1
  I.ME[2,3,3] <- 1
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  # per-SNP probability of mismatch
  P.Mismatch <- matrix(NA, 3, nSnp,
                       dimnames = list(c("DUP", "OH", "ME"), seq_along(MAF)))
  # probably possible to do it slightly faster w/o for loop, but not easily.
  for (x in seq_along(MAF)) {
    q <- MAF[x]
    AHWE <- c((1-q)^2, 2*q*(1-q), q^2)

    # duplicate mismatch
    P.Mismatch["DUP", x] <- 1 - sum(ErrA[x,,]^2 %*% AHWE)

    # parent-offspring opposite homozygosite
    AKAP <- matrix(c(1-q, (1-q)/2, 0,
                     q,   1/2,     1-q,
                     0,   q/2,     q), nrow=3, byrow=TRUE)
    P.AA <- sweep(AKAP, 2, AHWE, "*")   # joined prob actual genotypes
    P.OO <- t(ErrA[x,,]) %*% P.AA %*% ErrA[x,,]     # joined prob obs genotypes
    P.Mismatch["OH", x] <- P.OO[1,3] + P.OO[3,1]

    # offspring-dam-sire mendelian error
    P.AAA <- sweep( sweep(AKA2P, 2, AHWE, "*"), 3, AHWE, "*")
    P.OOA <- array(dim=c(3,3,3))
    for (i in 1:3) {
      P.OOA[,,i] <- t(ErrA[x,,]) %*% P.AAA[,,i]  %*% ErrA[x,,]
    }
    P.OOO <- array(dim=c(3,3,3))  # 3x obs
    for (j in 1:3) {
      P.OOO[,j,] <- P.OOA[,j,] %*% ErrA[x,,]
    }
    P.Mismatch["ME", x] <- sum(P.OOO * I.ME)
  }

  # rounding (especially necessary if Err=0.0)
  max.digits <- floor(abs(log10(.Machine$double.ep)))
  P.Mismatch <- round(P.Mismatch, digits = max.digits)

  if (Return == 'Counts') {
    MaxMismatchV <- setNames(numeric(3), c("DUP", "OH", "ME"))
    for (x in 1:3) {
      MaxMismatchV[x] <- stats::qbinom(qntl, size = nSnp, prob = mean(P.Mismatch[x, ]))
    }

    return( MaxMismatchV )
  } else {
    return(setNames(P.Mismatch, c("DUP", "OH", "ME")))
  }
}

#===============================================================================
#===============================================================================
