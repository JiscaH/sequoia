#' @title Estimate genotyping error rate (REMOVED; will be re-implemented)
#'
#' @description Estimate the genotyping error rates in SNP data, based on a
#'   pedigree and/or duplicates. Estimates probabilities (observed given
#'   actual) hom|other hom, het|hom, and hom|het. THESE ARE APPROXIMATE VALUES!
#'
#' @param GenoM  Genotype matrix
#' @param Pedigree  data.frame with columns id - dam - sire
#' @param Duplicates  matrix or data.frame with 2 columns, id1 & id2
#' @param Er_start  vector of length 3 with starting values for \code{optim}.
#' @param perSNP  logical, estimate error rate per SNP. WARNING not very
#'   precise, use only as an approximate indicator! Try on simulated data first,
#'   e.g. with \code{\link{SimGeno}}.
#'
#' @return  vector of length 3 with estimated genotyping error rates: the
#'  probabilities that
#' \itemize{
#'  \item hom|hom: an actual homozygote is observed as the other homozygote
#'  \item het|hom: an actual homozygote is observed as heterozygote
#'  \item hom|het: an actual heterozygote is observed as homozygote
#'  }
#'
#'  These are three independent parameters, that define the genotyping error
#'  matrix (see \code{\link{ErrToM}}) as follows:
#'
#' \tabular{lccc}{
#'             \tab \strong{0}       \tab \strong{1}   \tab \strong{2}  \cr
#'  \strong{0} \tab \eqn{1-E_1-E_2} \tab \eqn{E_2}    \tab \eqn{E_1} \cr
#'  \strong{1} \tab \eqn{E_3}        \tab \eqn{1-2E_3} \tab \eqn{E_3}   \cr
#'  \strong{2} \tab \eqn{E_1}       \tab \eqn{E_2}    \tab \eqn{1-E_1-E_2} \cr
#' }
#'
#'  Note that for \code{optim} a lower bound of 1e-6 and upper bound of 0.499
#'  are used; if these values are returned this should be interpreted as
#'  'inestimably small' and 'inestimably large', respectively. PLEASE DO NOT USE
#'  THESE VALUES AS INPUT IN SUBSEQUENT ANALYSIS BUT SUBSITUTE BY A SENSIBLE
#'  VALUE!!
#'
#' @details The result should be interpreted as approximate, ballpark estimates!
#' The estimated error rates from a pedigree will not be as accurate as from
#' duplicate samples. Errors in individuals without parents or offspring will
#' not be counted, and errors in individuals with only few offspring may not be
#' noted either. Deviation of genotype frequencies among founders from
#' Hardy-Weinberg equilibrium may wrongly be attributed to genotyping errors.
#' Last but not least, any pedigree errors will result in higher estimated
#' genotyping errors.
#'
#'
#' @importFrom stats optim setNames na.exclude
#'
#' @examples
#' GenoX <- SimGeno(Ped_griffin, nSnp=400, SnpError=c(0.01,0.07, 0.1),
#'                 ParMis=0.1, CallRate=0.9)
#' # EstEr(GenoM=GenoX, Pedigree=Ped_griffin)
#'
#' @keywords internal 
#' @noRd

EstEr <- function(GenoM, Pedigree, Duplicates=NULL, Er_start=c(.05,.05,.05),
                  perSNP=FALSE)
{

  stop("This function has been removed and will be re-implemented")

}

