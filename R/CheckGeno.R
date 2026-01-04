#=======================================================================
#' @title Check Genotype Matrix
#'
#' @description Check that the provided genotype matrix is in the correct
#'   format, and check for low call rate samples and SNPs.
#'
#' @section Thresholds: Appropriate call rate thresholds for SNPs and
#'   individuals depend on the total number of SNPs, distribution of call rates,
#'   genotyping errors, and the proportion of candidate parents that are SNPd
#'   (sibship clustering is more prone to false positives). Note that filtering
#'   first on SNP call rate tends to keep more individuals in.
#'
#' @param GenoM the genotype matrix.
#' @param quiet suppress messages.
#' @param Plot  display the plots of \code{\link{SnpStats}}.
#' @param Return  either 'GenoM' to return the cleaned-up genotype matrix, or
#'   'excl' to return a list with excluded SNPs and individuals (see Value).
#' @param DumPrefix length 2 vector, to check if these don't occur among
#'   genotyped individuals.
#' @param Strict Exclude any individuals genotyped for <5% of SNPs, and any SNPs
#'   genotyped for <5% of individuals (TRUE); this was the unavoidable default
#'   up to version 2.4.1. Otherwise only excluded are (very nearly) monomorphic
#'   SNPs, SNPs scored for fewer than 2 individuals, and individuals scored for
#'   fewer than 2 SNPs.
#'
#' @return If \code{Return='excl'} a list with, if any are found:
#'  \item{ExcludedSNPs}{SNPs scored for <10% of individuals; automatically
#'    excluded when running \code{\link{sequoia}}}
#'  \item{ExcludedSnps-mono}{monomorphic (fixed) SNPs; automatically excluded
#'    when running \code{\link{sequoia}}. This includes nearly-fixed SNPs with
#'    MAF \eqn{= 1/2N}. Column numbers are *after* removal of
#'    \code{ExcludedSNPs}, if any.}
#'  \item{ExcludedIndiv}{Individuals scored for <5% of SNPs; these cannot be
#'    reliably included during pedigree reconstruction. Individual call rate is
#'    calculated after removal of 'Excluded SNPs'}
#'  \item{Snps-LowCallRate}{SNPs scored for 10% -- 50% of individuals; strongly
#'    recommended to be filtered out}
#'  \item{Indiv-LowCallRate}{individuals scored for <50% of SNPs; strongly
#'    recommended to be filtered out}
#'
#' When \code{Return='excl'} the return is \code{\link{invisible}}, i.e. a check
#' is run and warnings or errors are always displayed, but nothing may be
#' returned.
#'
#' @seealso \code{\link{SnpStats}} to calculate SNP call rates;
#'   \code{\link{CalcOHLLR}} to count the number of SNPs scored in both focal
#'   individual and parent.
#'
#' @examples
#' GenoM <- SimGeno(Ped_HSg5, nSnp=400, CallRate = runif(400, 0.2, 0.8))
#' # the quick way:
#' GenoM.checked <- CheckGeno(GenoM, Return="GenoM")
#'
#' # the user supervised way:
#' Excl <- CheckGeno(GenoM, Return = "excl")
#' GenoM.orig <- GenoM   # make a 'backup' copy
#' if ("ExcludedSnps" %in% names(Excl))
#'   GenoM <- GenoM[, -Excl[["ExcludedSnps"]]]
#' if ("ExcludedSnps-mono" %in% names(Excl))
#'   GenoM <- GenoM[, -Excl[["ExcludedSnps-mono"]]]
#' if ("ExcludedIndiv" %in% names(Excl))
#'   GenoM <- GenoM[!rownames(GenoM) %in% Excl[["ExcludedIndiv"]], ]
#'
#' # warning about  SNPs scored for <50% of individuals ?
#' # note: this is not necessarily a problem, and sometimes unavoidable.
#' SnpCallRate <- apply(GenoM, MARGIN=2,
#'                      FUN = function(x) sum(x!=-9)) / nrow(GenoM)
#' hist(SnpCallRate, breaks=50, col="grey")
#' GenoM <- GenoM[, SnpCallRate > 0.6]
#'
#' # to filter out low call rate individuals: (also not necessarily a problem)
#' IndivCallRate <- apply(GenoM, MARGIN=1,
#'                        FUN = function(x) sum(x!=-9)) / ncol(GenoM)
#' hist(IndivCallRate, breaks=50, col="grey")
#' GoodSamples <- rownames(GenoM)[ IndivCallRate > 0.8]
#'
#' @export

CheckGeno <- function(GenoM, quiet=FALSE, Plot=FALSE,
                      Return="GenoM", Strict=TRUE, DumPrefix = c("F0", "M0"))
{

  if (!(isTRUE(quiet) | isFALSE(quiet)))  stop("'quiet' must be TRUE or FALSE")

  # basic checks ----
  if (is.null(GenoM)) stop("please provide 'GenoM'")
  if (!is.matrix(GenoM)) stop("'GenoM' should be a numeric matrix")
  GenoM[is.na(GenoM)] <- -9
  if (!all(GenoM %in% c(0,1,2,-1, -9))) {
    UniqueValues <- unique(c(GenoM))
    InvalidValues <- UniqueValues[!UniqueValues %in% c(0,1,2,-1,-9)]
    stop(paste0("'GenoM' includes invalid values: ", "'",
                paste(InvalidValues, collapse="', '"), "'"))
  }
  if (!Return %in% c("GenoM", "excl"))  stop("'Return' must be 'GenoM' or 'excl'")


  # check IDs ----
  if (is.null(rownames(GenoM)))
    stop("'GenoM' has no rownames, these should be the individual IDs")
  if (any(duplicated(rownames(GenoM))))
    stop("'GenoM' has duplicate IDs. Please exclude or rename these samples,",
         " or run GenoConvert() with UseFID=TRUE.")
  if (any(grepl(" ", rownames(GenoM))))
    stop("GenoM IDs (rownames) must not include spaces")

  DP <- sapply(DumPrefix, function(x) ifelse(nchar(x)==1, paste0(x,"0"), x))
  for (i in seq_along(DumPrefix)) {
    if (any(substr(rownames(GenoM),1,nchar(DP[i]))==DP[i])) {
      stop("GenoM IDs (rownames) may not start with Dummy prefix+0,",
           " use argument DummyPrefix to change them (current: ", DP[1], ', ', DP[2], ")")
    }
  }

  # Exclude low quality SNPs ----
  Excl <- list()
  sstats <- SnpStats(GenoM, Plot = Plot)

  if (Strict) {
    SNPs.TooMuchMissing <- sstats[,"Mis"] >= 0.95
  } else {
    SNPs.TooMuchMissing <- sstats[,"Mis"] >= 1 - 1/nrow(GenoM)  # genotyped for 0 or 1 indiv --> no use
  }
  if (any(SNPs.TooMuchMissing)) {
    if (!quiet)  cli::cli_alert_warning(c("There are {sum(SNPs.TooMuchMissing)}",
         " SNPs scored for {ifelse(Strict, '<5% of', '0 or 1')} individuals, these will be excluded"))
    GenoM <- GenoM[, !SNPs.TooMuchMissing]
    Excl[["ExcludedSnps"]] <- which(SNPs.TooMuchMissing)
    sstats <- SnpStats(GenoM, Plot = FALSE)   # update SnpStats
  }

  SNPs.mono <- sstats[,"AF"] <= 1/(2*nrow(GenoM)) | sstats[,"AF"] >= 1 - 1/(2*nrow(GenoM))
  if (any(SNPs.mono)) {
    if (!quiet)  cli::cli_alert_warning("There are {sum(SNPs.mono)} monomorphic (fixed) SNPs, these will be excluded")
    GenoM <- GenoM[, !SNPs.mono]
    Excl[["ExcludedSnps-mono"]] <- which(SNPs.mono)
  }

  SNPs.LotMissing <- apply(GenoM, 2, function(x) sum(x>=0)) < nrow(GenoM)/2
  if (any(SNPs.LotMissing)) {
    if (!quiet)  cli::cli_alert_warning(c('{ifelse("ExcludedSnps" %in% names(Excl), "In addition, there", "There")}',
         " are {sum(SNPs.LotMissing)} SNPs scored for <50% of individuals"))
    Excl[["Snps-LowCallRate"]] <- which(SNPs.LotMissing)
  }


  # Exclude low quality individuals ----
  Lscored <- apply(GenoM, 1, function(x) sum(x>=0))
  if (Strict) {
    Indiv.TooMuchMissing <- Lscored < ncol(GenoM)/20
  } else {
    Indiv.TooMuchMissing <- Lscored <= 1  # genotyped for 0 or 1 SNPs --> no point
  }
  if (any(Indiv.TooMuchMissing)) {
    cli::cli_alert_danger(c("There are {.strong {sum(Indiv.TooMuchMissing)} individuals}",
         ' scored for {ifelse(Strict, "<5% of", "0 or 1")} SNPs, ',
         "these {.strong WILL BE IGNORED}"))
    Excl[["ExcludedIndiv"]] <- rownames(GenoM)[which(Indiv.TooMuchMissing)]
    GenoM <- GenoM[!Indiv.TooMuchMissing, ]
  }

  Indiv.LotMissing <- apply(GenoM, 1, function(x) sum(x>=0)) < ncol(GenoM)/5
  if (any(Indiv.LotMissing)) {
     cli::cli_alert_danger(c('{ifelse("ExcludedIndiv" %in% names(Excl), "In addition, there", "There")}',
         " are {sum(Indiv.LotMissing)} individuals scored for <20% of SNPs,",
         "it is advised to treat their assignments with caution"))
    Excl[["Indiv-LowCallRate"]] <- rownames(GenoM)[Indiv.LotMissing]
  }


  if (!quiet) {
    MSG <- paste("There are ", nrow(GenoM), " individuals and ", ncol(GenoM), " SNPs.\n")
    if (any(c("ExcludedSnps", "ExcludedSnps-mono", "ExcludedIndiv") %in% names(Excl))) {
      cli::cli_alert_info(c("After exclusion, ", MSG))
    } else if (length(Excl)>0) {
      cli::cli_alert_info(MSG)
    } else {
      cli::cli_alert_success(c("Genotype matrix looks OK! ", MSG))
    }
  }


  if (Return == "GenoM") {
    return( GenoM )
  } else {
    invisible( Excl )
  }
}

