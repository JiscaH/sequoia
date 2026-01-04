#=======================================================================
#' @title Calculate assignment probabilities
#'
#' @description For each assigned offspring-parent pair, calculate the
#'  probability they are parent-offspring vs otherwise related. Probabilities
#'  are scaled to sum to one across all possible* relationships between the pair
#'  or trio; see Details.
#'
#' @details The returned probabilities are calculated from the likelihoods used
#'   throughout the rest of this package, by scaling them to sum to one across
#'   all possible relationships. For \code{Complex='simp'} these are
#'   PO=parent-offspring, FS=full siblings, HS=half siblings, GP=grand-parental,
#'   FA=full avuncular, HA=third degree relatives (incl half avuncular), and
#'   U=unrelated. For \code{Complex='full'} there are numerous double
#'   relationship considered (PO & HS, HS & HA, etc), making both numerator and
#'   denominator in the scaling step less unambiguous, and the returned
#'   probabilities an approximation.
#'
#'  The likelihoods are calculated by calling \code{\link{CalcPairLL}} once or
#'  twice for each id-dam and id-sire pair: once not conditioning on the
#'  co-parent, and once conditional on the co-parent, if any. For genotyped
#'  individuals this is done with \code{focal='PO'}, and for dummy individuals
#'  with \code{focal='GP'}.
#'
#'  For relationships between a genotyped and a dummy individual, it may only be
#'  possible to determine that the genotyped individual is a second degree
#'  relative (GP, HS, or FA) to the dummy's offspring. This then results in a
#'  probability of at most 0.33, even when the two are indeed parent and
#'  offspring.
#'
#'  See \code{\link{CalcPairLL}} and the vignettes for further details.
#'
#' Note that for large pedigrees this function can be fairly slow, especially
#' when using \code{\link{CalcPairLL}}'s default \code{Module='ped'} and
#' \code{Complex='full'}.
#'
#' Subsetting the genotype data may give different results, as the likelihoods
#'  and thus the probabilities depend on the allele frequencies in the sample.
#'
#' @param Pedigree  dataframe with columns id-dam-sire. By default, any
#'   non-genotyped individuals are 'dummified'; use \code{Module='par'} to
#'   ignore them.
#' @param GenoM  numeric matrix with genotype data: One row per individual,
#'   one column per SNP, coded as 0, 1, 2, missing values as a negative number
#'   or NA. You can reformat data with \code{\link{GenoConvert}}, or use other
#'   packages to get it into a genlight object and then use \code{as.matrix}.
#' @param nCores number of computer cores to use. If \code{2} or \code{4},
#'   package \pkg{parallel} is used (other values are not applicable).
#' @param quiet logical, suppress messages. No progress is printed when >1 core
#'   is used.
#' @param ... Additional arguments passed to \code{\link{CalcPairLL}}, such as
#'   the genotyping error rate \code{Err}, age information in
#'   \code{LifeHistData} and \code{AgePrior}, or \code{InclDup} to include the
#'   probability that the two samples are duplicates.
#'
#' @return the \code{Pedigree} dataframe with the three applicable columns
#' renamed to id-dam-sire, and 7 additional columns:
#' \item{Probdam}{Probability that individual in dam column is the maternal
#'   parent, rather than otherwise related (LL(PO)/sum(LL))}
#' \item{Probsire}{Analogous for sire}
#' \item{Probpair}{Probability for id-dam-sire trio. Approximated as the minimum
#'   of dam conditional on sire and sire conditional on dam, thus not including
#'   e.g. both being siblings (those other configurations are considered by
#'   sequoia during pedigree reconstruction, but can (currently) not be accessed
#'   directly)}
#' \item{dam_alt, sire_alt}{Most likely alternative (not PO) relationship
#'  between id-dam and id-sire, respectively}
#' \item{Probdam_alt, Probsire_alt}{Probability of most likely alternative
#'  relationship}
#'
#' @section Warning:
#' The probabilities will be less reliable with close inbreeding and double
#' relationships. This function has not been tested yet with hermaphrodites, and
#' is unlikely to give reliable results without further code updates.
#'
#' @seealso  \code{\link{CalcPairLL}}, \code{\link{LLtoProb}}
#'
#' @examples
#' test_ped <- Ped_griffin[21:25,]
#' # add an incorrect sire to illustrate
#' test_ped$sire <- as.character(test_ped$sire)
#' test_ped$sire[5] <- 'i057_2003_M'
#' Ped_with_probs <- CalcParentProbs(test_ped, Geno_griffin)
#' print(Ped_with_probs, digits=2)
#' # Any non-genotyped non-'dummifiable' individuals are automatically skipped
#'
#' # To get likelihoods for 'all' relationships, not just probabilities for
#' # PO & (next-)most-likely:
#' LL_sire_single <- CalcPairLL(
#'   Pairs = data.frame(id1=test_ped$id,
#'                      id2=test_ped$sire,
#'                      dropPar1='both', # drop both -> id2 as single parent
#'                      focal='PO'),
#'   Pedigree = Ped_griffin,   # pedigree to condition on
#'   GenoM = Geno_griffin, Plot=FALSE)
#' @export

CalcParentProbs <- function(Pedigree = NULL, GenoM = NULL, quiet = FALSE, nCores = 1, ...)

{
  #=========================
  # input checks ----
  if (!(isTRUE(quiet) | isFALSE(quiet)))  stop("'quiet' must be TRUE or FALSE")

  ## check genotype data ----
  GenoM <- CheckGeno(GenoM, quiet=TRUE, Plot=FALSE)
  # Do NOT remove any IDs from GenoM that are not in pedigree:
  # this affects the allele frequencies, which affects the likelihoods.

  ## check pedigree ---
  # drop individuals not in GenoM from Pedigree
  gID <- rownames(GenoM)
  # if Module='par', remove all non-SNPd individuals from the pedigree
  if (is.null(match.call()[['Module']]))  Module <- 'ped'  # default for CalcPairLL()
  Ped <- PedPolish(Pedigree, gID=gID, DropNonSNPd=(Module=='par'),
                   NullOK=FALSE, KeepAllColumns=FALSE)
  utils::flush.console()

  if (Module=='ped') keep_cat <- c('G','D') else keep_cat <- 'G'
  PedC <- getAssignCat(Ped, gID, minSibSize = "1sib1GP")
  PedC <- PedC[PedC$id.cat %in% keep_cat, ]


  #=========================
  # prepair pairs ----
  PairsL_IN <- list(dam_solo = with(PedC[PedC$dam.cat %in% keep_cat,],
                                    data.frame(ID1 = id, ID2 = dam, Sex2 = 1,
                                               patmat = 1, dropPar1 = 'both')),  # focus: PO, via maternal side
                    sire_solo = with(PedC[PedC$sire.cat %in% keep_cat,],
                                     data.frame(ID1 = id, ID2 = sire, Sex2 = 2,
                                                patmat = 2, dropPar1 = 'both')),  # focus: PO, via paternal side
                    dam_paired = with(PedC[PedC$dam.cat %in% keep_cat & PedC$sire.cat %in% keep_cat,],
                                      data.frame(ID1 = id, ID2 = dam, Sex2 = 1,
                                                 patmat = 1, dropPar1 = 'dam')),
                    sire_paired = with(PedC[PedC$dam.cat %in% keep_cat & PedC$sire.cat %in% keep_cat,],
                                       data.frame(ID1 = id, ID2 = sire, Sex2 = 2,
                                                  patmat = 2, dropPar1 = 'sire'))
  )
  for (x in names(PairsL_IN)) {
    PairsL_IN[[x]] <- cbind(PairsL_IN[[x]],
                            Sex1 = NA, dropPar2 = 'none',
                            focal = ifelse(PairsL_IN[[x]]$ID1 %in% gID, 'PO', 'GP'))
    # if ID1 is dummy, rel is relative to sibship.
    # Use U rather than GP to also calc FS/HS between dummies with assignment of
    # ID2's parent to ID1 as side-effect. These are excluded during GP assignment (focal=GP).
  }


  #=========================
  # call CalcPairLL ----
  # CalcPairLL() calls Fortran to do the calculations
  args_pairLL <- list(...)   # all optional arguments
  if (!nCores %in% c(1,2,4))  stop("nCores must be 1, 2, or 4")
  if (nCores == 1) {
    PairsL_LL <- list()
    for (x in names(PairsL_IN)) {
      if (!quiet) message(' ')
      if (!quiet) cli::cli_alert_info(paste("Calculating LL for ", x))
      PairsL_LL[[x]] <- CalcPairLL(Pairs = PairsL_IN[[x]],
                                   GenoM = GenoM, Pedigree = Ped, quiet=quiet,
                                   Plot=FALSE, ...)
    }

  } else {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      if (interactive() & !quiet) {
        cli::cli_alert_info("Installing package {.pkg parallel} to run on multiple cores... ")
      }
      utils::install.packages("parallel")
    }
    utils::flush.console()

    ARGS <- c(list(GenoM = GenoM, Pedigree = Ped, quiet=TRUE, Plot=FALSE),
              args_pairLL)

    cl <- parallel::makeCluster(nCores)
    PairsL_LL <- parallel::parLapply(cl, X=PairsL_IN,
                                     fun = function(PairsIN, ARGS) {
                                       do.call(sequoia::CalcPairLL,
                                               c(list(Pairs=PairsIN), ARGS)) },
                                     ARGS)

    parallel::stopCluster(cl)
    names(PairsL_LL) <- names(PairsL_IN)
  }


  #=========================
  # reorder for dummies ----
  # see CalcPairLL() help
  RelNames <- c("PO", "FS", "HS", "GP", "FA", "HA", "U")
  if (isTRUE(args_pairLL[['InclDup']]))  RelNames <- c('DUP', RelNames)

  if (is.null(match.call()[['Complex']])) {
    cmplx <- 'full'  # default for CalcPairLL()
  } else {
    cmplx <- match.call()[['Complex']]  # 'Complex' matches to base R function
  }

  for (x in names(PairsL_IN)) {
    isDum <- !PairsL_LL[[x]]$ID1 %in% gID   # ID1 is dummy
    if (any(isDum)) {
      PairsL_LL[[x]][isDum, RelNames] <- ReOrderDums( PairsL_LL[[x]][isDum, RelNames],
      Complex=cmplx)
    }
  }


  #=========================
  # scale likelihoods to probabilities ----
  PairsL_probs <- list()
  for (x in names(PairsL_IN)) {
    PairsL_probs[[x]] <- cbind(PairsL_LL[[x]][, c('ID1','ID2')],
                               plyr::aaply(as.matrix(PairsL_LL[[x]][,RelNames]),
                                           .margin=1, .fun=LLtoProb))
  }


  #=========================
  # identify second most likely relationship ----
  # not for pair: only subset of possible trio relationships is considered here.
  RelNames_noPO <- setdiff(RelNames, 'PO')
  for (x in c('dam_solo', 'sire_solo')) {
    M <- as.matrix(PairsL_probs[[x]][,  RelNames_noPO])
    wmx <- plyr::aaply(M, .margins=1, .fun=function(x) ifelse(all(is.na(x)), NA, which.max(x)))
    PairsL_probs[[x]]$rel_alt <-  RelNames_noPO[ wmx ]
    PairsL_probs[[x]]$prob_alt <- M[ cbind(1:nrow(M), wmx) ]
  }


  #=========================
  # add probabilities to pedigree ----
  cols_select <- list(solo = c('ID1', 'ID2', 'PO', 'rel_alt', 'prob_alt'),
                      paired = c('ID1', 'ID2', 'PO'))
  newnames <- list(dam = c('id', 'dam', 'Probdam', 'dam_alt', 'Probdam_alt', 'Probpair_dam'),
                   sire = c('id', 'sire', 'Probsire', 'sire_alt', 'Probsire_alt', 'Probpair_sire'))

  Probs_dam <- merge(PairsL_probs[['dam_solo']][, cols_select[['solo']]],
                     PairsL_probs[['dam_paired']][, cols_select[['paired']]],
                     by=c('ID1','ID2'), all.x=TRUE)
  names(Probs_dam) <-newnames[['dam']]
  Probs_sire <- merge(PairsL_probs[['sire_solo']][, cols_select[['solo']]],
                     PairsL_probs[['sire_paired']][, cols_select[['paired']]],
                     by=c('ID1','ID2'), all.x=TRUE)
  names(Probs_sire) <- newnames[['sire']]

  ped_colnames <- colnames(Pedigree)
  Pedigree$rownum <- 1:nrow(Pedigree)  # merge does not respect sort=FALSE
  Ped_out <- merge(Pedigree, Probs_dam, all.x=TRUE)
  Ped_out <- merge(Ped_out, Probs_sire, all.x=TRUE)
  # fix row & column order
  Ped_out <- Ped_out[Ped_out$rownum,
                     c(ped_colnames, newnames[['dam']][-c(1:2)], newnames[['sire']][-c(1:2)])]

  # for pair, minimum of dam-with-sire and sire-with-dam
  # (identical to LLRpair in sequoia() )
  Ped_out$Probpair <- with(Ped_out, pmin(Probpair_dam, Probpair_sire, na.rm=TRUE))
  Ped_out$Probpair_dam <- NULL  # drop column
  Ped_out$Probpair_sire <- NULL

  return(Ped_out)
}



#=======================================================================
################################################################################
#=======================================================================

#' @title transform log-likelihoods to probabilities
#'
#' @description transform a vector with log10 likelihoods to a vector with
#'  probabilities summing to one.
#'
#' @details The returned probabilities are calculated from the likelihoods used
#'   throughout the rest of this package, by scaling them to sum to one across
#'   all possible relationships. For \code{Complex='simp'} these are
#'   PO=parent-offspring, FS=full siblings, HS=half siblings, GP=grand-parental,
#'   FA=full avuncular, HA=third degree relatives (incl half avuncular), and
#' U=unrelated. For \code{Complex='full'} there are numerous double relationship
#'   considered (PO & HS, HS & HA, etc), making both numerator and denominator
#'   in the scaling step less unambiguous, and the returned probabilities an
#'   approximation.
#'
#'  Computational under/overflow issues are reduced by subtracted the
#'   maximum value before converting from log to regular scale. Probabilities
#'   that would still be smaller than the machine precision (\code{(LL -
#'   min(LL)/2) < log10(.Machine$double.xmin)}) are set to NA en then to 0,
#'   instead of \code{-Inf}, to avoid issues when scaling to sum to 1.
#'
#' @param LLv a vector with log10-likelihoods. All values >0 are set to NA.
#'
#' @return a vector with probabilities, with the same length and names.
#'
#' @examples
#' LL_pairs <- CalcPairLL(data.frame(ID1='i042_2003_F',
#'                          ID2=c('i015_2001_F', 'i022_2002_F', 'i035_2002_F')),
#'                   GenoM = Geno_griffin, Complex='simp', Err=1e-3, Plot=FALSE)
#' prob_pairs <- t(apply(LL_pairs[,10:16], MARGIN=1, LLtoProb))
#' # - or -
#' prob_pairs <- plyr::aaply(as.matrix(LL_pairs[,10:16]), .margin=1, LLtoProb)
#' round(prob_pairs, 3)
#'
#' # i035_2002_F is MHS of i042_2003_F, but when not conditioning on any other
#' # relatives has a higher LL to be 3rd degree relative (HA)
#' # (possibly genotyping errors, or just randomness of Mendelian inheritance)
#'
#' @export

LLtoProb <- function(LLv)   # vector with likelihoods
{
  # set various NA codes to NA
  LLv[LLv >0] <- NA
  # all LL's are equal in case of dummy parent with no GPs or further genotyped offspring
  # set to NA: no information.
  if (suppressWarnings(max(LLv, na.rm=TRUE) - min(LLv, na.rm=TRUE)) < 0.01) LLv[] <- NA
  # subtract the maximum LL from all
  # (i.e. divide on non-logscale, which does not change the ratios)
  LLv <- LLv - suppressWarnings(max(LLv, na.rm=TRUE))
  # check if can be converted to non-log without under/overflow issues
  machine_min <- log10(.Machine$double.xmin)
  min_dLL <- suppressWarnings(min(LLv, na.rm=TRUE))
  toolow <- NULL
  if (min_dLL < machine_min) {
    # set the smallest values to have a probability of 0 (rather than -Inf)
    toolow <- which(LLv < machine_min)
    LLv[toolow] <- NA
    LLv <- LLv - suppressWarnings(min(LLv, na.rm=TRUE))/2
  }
  probv <- 10**LLv
  probv <- probv / sum(probv, na.rm=TRUE)
  if (!is.null(toolow))  probv[toolow] <- 0.0

  return( probv )
}


#=======================================================================
################################################################################
#=======================================================================

#' @title Re-order likelihood vector
#'
#' @description change the order for dummy individuals so that the
#'  results are easily comparable with genotyped individuals
#'
#' @param LLM a matrix with log10-likelihoods; individuals in rows.
#'
#' @return a matrix with similar dimensions, but changed column order:
#'   \item{.}{GP \eqn{\rightarrow} PO (GP of sibship = parent of dummy)}
#'   \item{.}{FA \eqn{\rightarrow} FS}
#'   \item{.}{HA \eqn{\rightarrow} HS}
#'   \item{.}{FS \eqn{\rightarrow} FA (FS of sibs=offspring of dummy, NOT parent)}
#'   \item{.}{HS \eqn{\rightarrow} HA}
#'   \item{.}{U \eqn{\rightarrow} U}
#'   \item{.}{PO \eqn{\rightarrow} DUP (only if 'DUP' already among column names)}
#'
#' @keywords internal  
#' @noRd

ReOrderDums <- function(LLM, Complex='full')
{
  RelNames <- c("PO", "FS", "HS", "GP", "FA", "HA", "U", 'DUP')
  if (!all(colnames(LLM) %in% RelNames)) {
    cli::cli_alert_danger('Valid column names are {RelNames}')
    stop()
  }

  LLOUT <- LLM   # same dims & dimnames
  if ('DUP' %in% colnames(LLM))  LLOUT[,'DUP'] <- LLM[,'PO']
  LLOUT[,'PO'] <- LLM[,'GP']
  LLOUT[,'FS'] <- LLM[,'FA']
  if (Complex!='mono')  LLOUT[,'HS'] <- LLM[,'HA']
  LLOUT[,'GP'] <- NA   # calculated, but not returned separately (included under HA)
  LLOUT[,'FA'] <- LLM[,'FS']
  if (Complex!='mono')  LLOUT[,'HA'] <- LLM[,'HS']
  # U remains U.
  return(LLOUT)
}
