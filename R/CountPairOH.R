#' @title Count opposing homozygous SNPs between pairs of individuals
#'
#' @description Quick identification of likely parents, as the number of
#'   opposing homozygous (OH) SNPs is expected to be zero for parent- offspring
#'   pairs in absence of genotyping errors, and greater than zero for all other
#'   pairs.
#'
#' @param x  Either a matrix, dataframe or similar where the first two columns
#'   are individual IDs, or a vector with IDs. In the second case, you may
#'   provide \code{ID2}, and the output will be an ID1 x ID2 matrix; else the
#'   output will be an ID1 x ID1 matrix. Non-genotyped individuals are included
#'   in the results with all \code{NA}'s, and a warning.
#' @param ID2  optional second vector with IDs
#' @param GenoM  numeric matrix with genotype data: One row per individual, one
#'   column per SNP, coded as 0, 1, 2, missing values as a negative number or
#'   NA. Row names must be individual IDs, column names are ignored. You can
#'   reformat data with \code{\link{GenoConvert}}, or use other packages to get
#'   it into a genlight object and then use \code{as.matrix}.
#' @param max_OH  stop counting OH's for a pair if this value is reached, to
#'   reduce computation time. Ignored if negative value or equal to total number
#'   of SNPs.
#' @param quiet  suppress messages
#'
#' @details  Counting the number of opposing homozygous (OH) SNPs is much faster
#' than calculating likelihoods, and does not rely on an estimated genotyping
#' error rate. It can therefore be useful during quality control, or to help
#' figure out problems when assignment rate with \code{\link{sequoia}} is
#'  lower than expected.
#'
#' @seealso \code{\link{CalcPairLL}} to calculate likelihoods for pairs,
#'  \code{\link{CalcOHLLR}} to calculate OH for a pedigree,
#'  \code{\link{CalcMaxMismatch}} for calculation of the maximum OH used by
#'  \code{\link{sequoia}} to filter potential parent-offspring pairs.
#'
#' @examples
#' offspring_ids <- with(LH_HSg5, ID[BirthYear==2001])
#' candidate_father_ids <- with(LH_HSg5, ID[BirthYear==2000 & Sex==2])
#'
#' OH_matrix <- CountOH(offspring_ids, candidate_father_ids, GenoM=Geno_HSg5)
#'
#' hist(c(OH_matrix), breaks=c(0:40)-.5)
#' # with high quality SNP data, there is often a clear separation in OH counts
#' # between parent-offspring pairs (here: OH<3) and others (here: OH>4).
#' # BUT: non-PO close relatives may have very low OH counts by chance,
#' # and true PO pairs may have fairly high OH counts due to genotyping errors.
#'
#'
#' @useDynLib sequoia, .registration = TRUE
#'
#' @export

CountOH <- function(x = NULL,
                    ID2 = NULL,
                    GenoM = NULL,
                    max_OH = -1,
                    quiet = FALSE)
{
  if (!(isTRUE(quiet) | isFALSE(quiet)))  stop("'quiet' must be TRUE or FALSE")
  if (length(dim(x)) == 2) {  # x is matrix or dataframe  or similar
    OutFormat = 'Dataframe'
    ID1 <- as.character(x[[1]])
    ID2 <- as.character(x[[2]])
  } else if (is.null(x)) {  # TODO: check that x is a vector
    stop('Please provide either a matrix/dataframe or two vectors with IDs')
  } else {
    OutFormat <- 'Matrix'
    ID1 <- as.character(x)
    if (is.null(ID2)) {
      if (!quiet)  cli::cli_alert_info('No ID2 provided, so using ID2 = ID1')
      ID2 <- as.character(x)
    } else {
      ID2 <- as.character(ID2)
    }
  }
  if (!all(c(ID1,ID2) %in% rownames(GenoM)))
    if (!quiet)  cli::cli_alert_info('Some IDs are not among rownames of GenoM')

  GenoM <- CheckGeno(GenoM, quiet = TRUE, Plot = FALSE)

  if (OutFormat == 'Dataframe') {
    Pairs <- as.data.frame(x[,1:2])
    names(Pairs) <- c('ID1', 'ID2')
  } else {
    Pairs <- data.frame(ID1 = rep(ID1, each = length(ID2)),
                        ID2 = rep(ID2, times = length(ID1)))
  }

  gID <- rownames(GenoM)
  GenoNums <- setNames(seq_along(gID), gID)

  for (x in 1:2) {
    Pairs[, paste0('ID',x,'.num')] <- ifelse(Pairs[, paste0('ID',x)] %in% gID,
                                             GenoNums[Pairs[, paste0('ID',x)]], 0)

  }
  Np <- nrow(Pairs)
  if (max_OH < 0)  max_OH <- ncol(GenoM)

  TMP <- .Fortran(countpairoh,
                  ng = as.integer(nrow(GenoM)),
                  nm = as.integer(ncol(GenoM)),   # TODO: check for consistent naming
                  np = as.integer(Np),
                  maxoh = as.integer(max_OH),
                  genofr = as.integer(GenoM),
                  pairids = as.integer(c(Pairs$ID1.num, Pairs$ID2.num)),
                  ohrf = integer(Np))

  TMP$ohrf[ TMP$ohrf <0 ] <- NA   # either or both not genotyped

  if (OutFormat == 'Dataframe') {
    return(
      cbind(x, OH = TMP$ohrf )
    )

  } else {
    return(
      matrix(TMP$ohrf, byrow = TRUE,
             nrow=length(ID1), ncol=length(ID2),
             dimnames = list(ID1 = ID1, ID2 = ID2))
    )
  }

}
