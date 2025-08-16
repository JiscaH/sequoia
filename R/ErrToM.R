#' @title Generate Genotyping Error Matrix
#'
#' @description Make a vector or matrix specifying the genotyping error
#'   pattern, or a function to generate such a vector/matrix from a single
#'   value Err.
#'
#' with the probabilities of observed genotypes
#'   (columns) conditional on actual genotypes (rows), or return a function to
#'   generate such matrices (using a single value Err as input to that
#'  function).
#'
#' @details By default (\code{flavour} = "version2.9"), \code{Err} is
#'   interpreted as a locus-level error rate (rather than allele-level), and
#'   equals the probability that an actual heterozygote is observed as either
#'   homozygote (i.e., the probability that it is observed as AA = probability
#'   that observed as aa = \code{Err}/2). The probability that one homozygote is
#'   observed as the other is (\code{Err}/2\eqn{)^2}.
#'
#' The inbuilt 'flavours' correspond to the presumed and simulated error
#' structures, which have changed with sequoia versions. The most appropriate
#' error structure will depend on the genotyping platform; 'version0.9' and
#' 'version1.1' were inspired by SNP array genotyping while 'version1.3' and
#' 'version2.0' are intended to be more general.
#'
#' This function, and throughout the package, it is assumed that the two alleles
#' \eqn{A} and \eqn{a} are equivalent. Thus, using notation \eqn{P}(observed
#' genotype |actual genotype), that \eqn{P(AA|aa) = P(aa|AA)},
#' \eqn{P(aa|Aa)=P(AA|Aa)}, and \eqn{P(aA|aa)=P(aA|AA)}.
#'
#' \tabular{lccc}{
#'  \strong{version} \tab \strong{hom|hom} \tab \strong{het|hom} \tab
#' \strong{hom|het} \cr
#'  \strong{2.9}  \tab \eqn{(E/2)^2} \tab \eqn{E-(E/2)^2} \tab \eqn{E/2} \cr
#'  \strong{2.0}  \tab \eqn{(E/2)^2} \tab \eqn{E(1-E/2)}  \tab \eqn{E/2} \cr
#'  \strong{1.3}  \tab \eqn{(E/2)^2} \tab \eqn{E}   \tab \eqn{E/2} \cr
#'  \strong{1.1}  \tab \eqn{E/2}     \tab \eqn{E/2} \tab \eqn{E/2} \cr
#'  \strong{0.9}  \tab \eqn{0}       \tab \eqn{E}   \tab \eqn{E/2} \cr
#' }
#'
#' or in matrix form, Pr(observed genotype (columns) | actual genotype (rows)):
#'
#' \emph{version2.9:}
#' \tabular{lccc}{
#'     \tab \strong{0} \tab \strong{1} \tab \strong{2} \cr
#'  \strong{0}  \tab \eqn{1-E} \tab \eqn{E -(E/2)^2} \tab \eqn{(E/2)^2} \cr
#'  \strong{1}  \tab \eqn{E/2}       \tab \eqn{1-E}      \tab \eqn{E/2}    \cr
#'  \strong{2}  \tab \eqn{(E/2)^2}   \tab \eqn{E -(E/2)^2} \tab \eqn{1-E} \cr
#' }
#'
#' \emph{version2.0:}
#' \tabular{lccc}{
#'     \tab \strong{0} \tab \strong{1} \tab \strong{2} \cr
#'  \strong{0}  \tab \eqn{(1-E/2)^2} \tab \eqn{E(1-E/2)} \tab \eqn{(E/2)^2} \cr
#'  \strong{1}  \tab \eqn{E/2}       \tab \eqn{1-E}      \tab \eqn{E/2}    \cr
#'  \strong{2}  \tab \eqn{(E/2)^2}   \tab \eqn{E(1-E/2)} \tab \eqn{(1-E/2)^2} \cr
#' }
#'
#' \emph{version1.3}
#' \tabular{lccc}{
#'     \tab \strong{0} \tab \strong{1} \tab \strong{2} \cr
#'  \strong{0}  \tab \eqn{1-E-(E/2)^2} \tab \eqn{E} \tab \eqn{(E/2)^2} \cr
#'  \strong{1}  \tab \eqn{E/2}       \tab \eqn{1-E}      \tab \eqn{E/2}    \cr
#'  \strong{2}  \tab \eqn{(E/2)^2}   \tab \eqn{E} \tab \eqn{1-E-(E/2)^2} \cr
#' }
#'
#' \emph{version1.1}
#' \tabular{lccc}{
#'     \tab \strong{0} \tab \strong{1} \tab \strong{2} \cr
#'  \strong{0}  \tab \eqn{1-E} \tab \eqn{E/2} \tab \eqn{E/2} \cr
#'  \strong{1}  \tab \eqn{E/2}       \tab \eqn{1-E}      \tab \eqn{E/2}    \cr
#'  \strong{2}  \tab \eqn{E/2}   \tab \eqn{E/2} \tab \eqn{1-E} \cr
#' }
#'
#' \emph{version0.9} (not recommended)
#' \tabular{lccc}{
#'     \tab \strong{0} \tab \strong{1} \tab \strong{2} \cr
#'  \strong{0}  \tab \eqn{1-E} \tab \eqn{E} \tab \eqn{0} \cr
#'  \strong{1}  \tab \eqn{E/2}       \tab \eqn{1-E}      \tab \eqn{E/2}    \cr
#'  \strong{2}  \tab \eqn{0}   \tab \eqn{E} \tab \eqn{1-E} \cr
#' }
#'
#' When \code{Err} is a length 3 vector, or if \code{Return = 'vector'} these
#'  are the following probabilities:
#' \itemize{
#'  \item hom|hom: an actual homozygote is observed as the other homozygote
#'   (\eqn{E_1})
#'  \item het|hom: an actual homozygote is observed as heterozygote (\eqn{E_2})
#'  \item hom|het: an actual heterozygote is observed as homozygote (\eqn{E_3})
#'  }
#'
#'  and Pr(observed genotype (columns) | actual genotype (rows)) is then:
#' \tabular{lccc}{
#'             \tab \strong{0}       \tab \strong{1}   \tab \strong{2}  \cr
#'  \strong{0} \tab \eqn{1-E_1-E_2} \tab \eqn{E_2}    \tab \eqn{E_1} \cr
#'  \strong{1} \tab \eqn{E_3}        \tab \eqn{1-2E_3} \tab \eqn{E_3}   \cr
#'  \strong{2} \tab \eqn{E_1}       \tab \eqn{E_2}    \tab \eqn{1-E_1-E_2} \cr
#' }
#'
#'  When the SNPs are scored via sequencing (e.g. RADseq or DArTseq), the 3rd
#'  error rate (hom|het) is typically considerably higher than the other two,
#'  while for SNP arrays it tends to be similar to P(het|hom).
#'
#' @param Err estimated genotyping error rate, as a single number, or 3x3 or 4x4
#'   matrix, or length 3 vector. If a single number, an error model is used that
#'   aims to deal with scoring errors typical for SNP arrays. If a matrix, this
#'   should be the probability of observed genotype (columns) conditional on
#'   actual genotype (rows). Each row must therefore sum to 1. If
#'   \code{Return='function'}, this may be \code{NA}. If a vector, these are the
#'   probabilities (observed given actual) hom|other hom, het|hom, and hom|het.
#' @param flavour vector-generating or matrix-generating function, or one of
#'   'version2.9', 'version2.0', 'version1.3' (='SNPchip'), 'version1.1'
#'   (='version111'), referring to the sequoia version in which it was used as
#'   default. Only used if \code{Err} is a single number.
#' @param Return output, 'matrix' (default), 'vector', 'function'
#'   (matrix-generating), or 'v_function' (vector-generating)
#'
#' @return Depending on \code{Return}, either:
#'  \itemize{
#'    \item \code{'matrix'}: a 3x3 matrix, with probabilities of observed genotypes
#'    (columns) conditional on actual (rows)
#'    \item \code{'function'}: a function taking a single value \code{Err} as input, and
#'    generating a 3x3 matrix
#'    \item \code{'vector'}: a length 3 vector, with the probabilities (observed given
#'      actual) hom|other hom, het|hom, and hom|het.
#'  }
#'
#'
#' @examples
#' ErM <- ErrToM(Err = 0.05)
#' ErM
#' ErrToM(ErM, Return = 'vector')
#'
#'
#' # use error matrix from Whalen, Gorjanc & Hickey 2018
#' funE <- function(E) {
#'  matrix(c(1-E*3/4, E/2, E/4,
#'           E/4, 1-2*E/4, E/4,
#'           E/4, E/2, 1-E*3/4),
#'           3,3, byrow=TRUE)  }
#' ErrToM(Err = 0.05, flavour = funE)
#' # equivalent to:
#' ErrToM(Err = c(0.05/4, 0.05/2, 0.05/4))
#'
#' @export

ErrToM <- function(Err = NA,
                   flavour = "version2.9",
                   Return = "matrix")
{
  if (is.null(Err) || length(Err)==1 && is.na(Err) && Return %in% c("function", "v_function"))
      Err <- 0.1   # only used for testing
  if (!is.atomic(Err) || !length(Err) %in% c(1,3,9,16))  stop("'Err' must be a single number, length 3 vector, or 3x3 matrix")
  if (any(Err<0 | Err>1) || !is.double(Err)) stop("'Err' must be (a) number(s) between 0 and 1")
  ErrV <- NULL  # hom|hom, het|hom, hom|het
  ErrM <- NULL  # observed (columns) conditional on actual (rows)
  ErFuncV <- NULL  # function to generate a length-3 vector
  ErFunc <- NULL   # function to generate a 3x3 or 4x4 matrix
  ValidFlavours <- c(paste0('version', c('2.9', '2.0', '1.3', '1.1', '0.9', '0.7')), 'SNPchip')

  if (is.matrix(Err)) {

    if (nrow(Err)==3 & ncol(Err)==3) {
      ErrM <- Err
    } else if (nrow(Err)==4 & ncol(Err)==4) {
      ErrM <- shrinkEM(Err)
    } else {
      stop("Error matrix should be a 3x3 or 4x4 matrix")
    }

  } else if (length(Err)==3) {
    ErrV <- Err

  } else if (length(Err) > 1) {
    stop("Err can be a 3x3 or 4x4 matrix, a length 3 vector, or a single number")
  } else {

    if (is.function(flavour)) {

      Err_test <- flavour(Err)
      if (length(Err_test)==3) {
        ErFuncV <- flavour
        ErrV <- ErFuncV(Err)
      } else {
        ErFunc <- flavour
        if (!is.matrix(Err_test))  stop("ErFunc(E) should return a 3x3 or 4x4 matrix")
        if (!(all(dim(Err_test)==4) | all(dim(Err_test)==3)) )  stop("ErFunc(E) should return a 4x4 or 3x3 matrix")
        ErrM <- ErFunc(Err)
      }
      Err_test.B <- flavour(Err+0.1)
      if (all(Err_test.B == Err_test))  stop("'flavour' function is not a function of error rate Err")

    } else if (length(flavour)==1 && flavour %in% ValidFlavours) {
      # hom|hom, het|hom, hom|het
      ErFuncV <- switch(flavour,
                        "version2.9" = function(E) c((E/2)^2, E-(E/2)^2, E/2),
                        "version2.0" = function(E) c((E/2)^2, E*(1-E/2), E/2),
                        "version1.3" = function(E) c((E/2)^2, E, E/2),
                        "SNPchip"    = function(E) c((E/2)^2, E, E/2),
                        "version1.1" = function(E) c(E/2, E/2, E/2),
                        "version0.9" = function(E) c(0, E, E/2),
                        "version0.7" = function(E) c(0, E, E/2),
                        function(E) c(0,0,0))

      ErrV <- ErFuncV(Err)
    } else {
      stop("Unknown ErrFlavour, choose 'version2.9', 'version2.0', 'version1.3', 'version1.1', \n",
           "or specify vector or matrix(-generating function) via 'Err'")
    }
  }

  if (!is.null(ErrM) && (!is.double(ErrM) || any(ErrM<0 | ErrM>1)))
    stop("Error vector/matrix values must be between 0 and 1")
  if (!is.null(ErrV) && (!is.double(ErrV) || any(ErrV<0 | ErrV>1)))
     stop("Error vector/matrix values must be between 0 and 1")

  if (!is.null(ErrM) && !all(abs(rowSums(ErrM) - 1) < sqrt(.Machine$double.eps)))
    stop("Error matrix rows must sum to 1")

  if (Return == "matrix") {
    if (is.null(ErrM))  ErrM <- ErV2M(ErrV)
    #dimnames(ErrM) <- list(paste0("act-", 0:2), paste0("obs-", 0:2, "|act"))
    return( ErrM )

  } else if (Return == 'vector') {
    if (is.null(ErrV))  ErrV <- ErM2V(ErrM)
    return( ErrV )

  } else if (Return == 'v_function') {
    if (is.null(ErFuncV) & !is.null(ErFunc))
      stop('Cannot automatically make vector function from matrix function')
    if (!is.null(ErFuncV)) {
      return( ErFuncV )
    } else {
      stop("Don't know how to make error function from vector or matrix")
    }
  } else if (Return == "function") {
    if (is.null(ErFunc) & !is.null(ErFuncV))  ErFunc <- EfuncV2M(ErFuncV)
    if (!is.null(ErFunc)) {
      return( ErFunc )
    } else {
      stop("Don't know how to make error function from vector or matrix")
    }

  } else {
    stop("Unknown Return format")
  }
}

#===============================================================================
#===============================================================================

# 4x4 matrix (aa, aA, Aa, AA) to 3x3 matrix
shrinkEM <- function(EM4) {
  EM3 <- matrix(NA, 3,3)
  EM3[c(1,3), c(1,3)] <- EM4[c(1,4), c(1,4)]
  EM3[2, c(1,3)] <- EM4[2, c(1,4)]+EM4[3, c(1,4)]
  EM3[c(1,3), 2] <- EM4[c(1,4), 2] + EM4[c(1,4), 2]
  EM3[2, 2] <- sum(EM4[2:3, 2:3])
  return( EM3 )
}


# vector -> matrix Pr(observed genotype (columns) | actual genotype (rows)):
ErV2M <- function(ErV,
                  CallRate=NULL) {  # for use by SimGeno
  if (is.null(CallRate)) {
    ErrM <- matrix(NA, 3,3, dimnames = list(act=0:2, obs=0:2))
  } else {
    ErrM <- matrix(NA, 3,4, dimnames = list(act=0:2, obs=-1:2))
  }
  ErrM['0', c('0','1','2')] <- c(1-ErV[1]-ErV[2],  ErV[2],  ErV[1])
  ErrM['1', c('0','1','2')] <- c(ErV[3],  1-2*ErV[3],  ErV[3])
  ErrM['2', c('0','1','2')] <- c(ErV[1],  ErV[2],  1-ErV[1]-ErV[2])

  if (!is.null(CallRate)) {
    ErrM[,c('0','1','2')] <- ErrM[,c('0','1','2')] * CallRate
    ErrM[,'-1'] <- 1 - CallRate
  }
  return( ErrM )
}


# matrix --> vector: hom -> other hom, hom -> het, het -> hom
ErM2V <- function(ErrM) {
  ErV <- setNames(rep(NA,3), c('hom|hom', 'het|hom', 'hom|het'))
  ErV['hom|hom'] <- (ErrM[1,3] + ErrM[3,1])/2
  ErV['het|hom'] <- (ErrM[1,2] + ErrM[3,2])/2
  ErV['hom|het'] <- (ErrM[2,1] + ErrM[2,3])/2
  return( ErV )
}


EfuncV2M <- function(Vfun) {
  VP <- body(Vfun)[2:4]
  Mfun <- function(E) {}
  body(Mfun) <- bquote(matrix(c(1-.(VP[[1]])-.(VP[[2]]), .(VP[[2]]), .(VP[[1]]),
                     .(VP[[3]]), 1-2*.(VP[[3]]), .(VP[[3]]),
                     .(VP[[1]]), .(VP[[2]]), 1-.(VP[[1]])-.(VP[[2]])),
                   3,3,byrow=TRUE))
  return(Mfun)
}
