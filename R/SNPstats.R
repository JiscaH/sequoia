#' @title SNP Summary Statistics
#'
#' @description Estimate allele frequency (AF), missingness and Mendelian
#' errors per SNP.
#'
#' @details Calculation of these summary statistics can be done in PLINK, and
#'   SNPs with low minor allele frequency or high missingness should be filtered
#'   out prior to pedigree reconstruction. This function is provided as an aid
#'   to inspect the relationship between AF, missingness and genotyping error to
#'   find a suitable combination of SNP filtering thresholds to use.
#'
#'   For pedigree reconstruction, SNPs with zero or one copies of the alternate
#'   allele in the dataset (MAF \eqn{\le 1/2N}) are considered fixed, and
#'   excluded.
#'
#' @param GenoM  genotype matrix, in sequoia's format: 1 column per SNP, 1 row
#'   per individual, genotypes coded as 0/1/2/-9, and row names giving individual
#'   IDs.
#' @param Pedigree  dataframe with 3 columns: ID - parent1 - parent2.
#'   Additional columns and non-genotyped individuals are ignored. Used to count
#'   Mendelian errors per SNP and (poorly) estimate the error rate.
#' @param Duplicates  dataframe with pairs of duplicated samples
#' @param Plot  logical, show histograms of the results?
#' @param quiet  logical, suppress messages?
#' @param calc_HWE logical, calculate chi-square test for Hardy-Weinberg
#'   equilibrium? Can be relatively time consuming for large datasets.
#' @param ErrFlavour DEPRECATED AND IGNORED. Was used to estimate \code{Err.hat}
#'
#' @return A matrix with a number of rows equal to the number of SNPs
#'  (=number of columns of GenoM), and when no Pedigree is provided 2 columns:
#' \item{AF}{Allele frequency of the 'second allele' (the one for which the
#'   homozygote is coded 2)}
#' \item{Mis}{Proportion of missing calls}
#' \item{HWE.p}{p-value from chi-square test for Hardy-Weinberg equilibrium}
#'
#' When a Pedigree is provided, there are 8 additional columns:
#' \item{n.dam, n.sire, n.pair}{Number of dams, sires, parent-pairs successfully
#'   genotyped for the SNP}
#' \item{OHdam, OHsire}{Count of number of opposing homozygous cases}
#' \item{MEpair}{Count of Mendelian errors, includes opposing homozygous cases
#'   when only one parent is genotyped}
#' \item{n.dups, n.diff}{Number of duplicate pairs successfully genotyped for
#'   the SNP; number of differences. The latter does not count cases where one
#'   duplicate is not successfully genotyped at the SNP}
#'
#' @seealso  \code{\link{GenoConvert}} to convert from various data formats;
#'   \code{\link{CheckGeno}} to check the data is in valid format for sequoia
#'   and exclude monomorphic SNPs etc., \code{\link{CalcOHLLR}} to calculate OH
#'   & ME per individual.
#'
#' @examples
#' Genotypes <- SimGeno(Ped_HSg5, nSnp=100, CallRate = runif(100, 0.5, 0.8),
#'                      SnpError = 0.05)
#' SnpStats(Genotypes)   # only plots; data is returned invisibly
#' SNPstats <- SnpStats(Genotypes, Pedigree=Ped_HSg5)
#'
#' @importFrom graphics plot points legend
#'
#' @export

SnpStats <- function(GenoM,
                     Pedigree = NULL,
                     Duplicates = NULL,
                     Plot = TRUE,
                     quiet = TRUE,
                     calc_HWE = TRUE,
                     ErrFlavour)
{
  if (!inherits(GenoM, 'matrix'))  stop("GenoM must be a matrix")
  GenoM[is.na(GenoM)] <- -9
  # missingness
  Mis <- apply(GenoM, 2, function(x) sum(x < 0))/nrow(GenoM)
  # allele frequency
  AF <- apply(GenoM, 2, function(x) sum(x[x>=0])/(2*sum(x>=0)))
  OUT <- cbind(AF, Mis)

  # Hardy-weinberg equilibrium
  if (!calc_HWE %in% c(TRUE, FALSE))  stop('calc_HWE must be TRUE or FALSE')
  if (calc_HWE) {
    counts.o <- apply(GenoM, 2, function(x) table(factor(x, levels=c(0,1,2))))  # observed genotype counts
    freq.e <- rbind('0' = (1-AF)^2,
                    '1' = 2*AF*(1-AF),
                    '2' = AF^2)
    freq.e[, Mis==1] <- 0.0
    calc.HWE.p <- function(i) {
      if (all(counts.o[,i]==0)) {
        NA  # all missing
      } else {
        suppressWarnings(stats::chisq.test(x = counts.o[,i],
                                    p = freq.e[,i])$p.value)
      }
    }
    HWE.p <- sapply(1:ncol(GenoM), calc.HWE.p)

    OUT <- cbind(OUT, HWE.p)
  }

  if (is.logical(Pedigree)) {
    Plot <- Pedigree
    Pedigree <- NULL
  }

  if (!is.null(Pedigree)) {
    Par <- PedPolish(Pedigree, gID = rownames(GenoM), DropNonSNPd=TRUE,
                     KeepAllColumns=FALSE)
    Par <- Par[!(is.na(Par[,'dam']) & is.na(Par[,'sire'])), ]
    if (nrow(Par)==0) {
      cli::cli_alert_warning('Cannot count OH, because pedigree does not have any genotyped parent-offspring pairs')
    } else if (!quiet) {
      N <- c(i=nrow(Par), d=sum(Par$dam %in% rownames(GenoM)), s=sum(Par$sire %in% rownames(GenoM)))
      cli::cli_alert_info("Conditioning on genotyped-only pedigree with {N['i']} individuals, {N['d']} dams and {N['s']} sires")
    }
  } else {
    if (!quiet)  cli::cli_alert_info("Not conditioning on any pedigree")
  }

  if (!is.null(Pedigree) | !is.null(Duplicates)) {
    OHcounts <- OHperSNP(GenoM, Par, Duplicates)
    OUT <- cbind(OUT, OHcounts)
  }

  if (Plot) {
    img <- tryCatch(
      {
        suppressWarnings( PlotSnpStats( OUT ) )
      },
      error = function(e) {
        message("Plotting area too small for SnpStats() plot (or other plotting problem)")
      })
  }

  if (is.null(colnames(GenoM)) | all(colnames(GenoM) == paste0('V', seq_len(ncol(GenoM))))) {
    w <- nchar(ncol(GenoM))
    rownames(OUT) <- paste0("SNP", formatC(1:nrow(OUT), width=w, flag="0"))
  } else {
    rownames(OUT) <- colnames(GenoM)  # SNP names
  }
  invisible( OUT )
}



#==============================================================================
#' @title Estimate Genotyping Error Rate
#'
#' @description Estimate genotyping error rate from Mendelian errors per SNP.
#'
#' @param GenoM  genotype matrix, in sequoia's format: 1 column per SNP, 1 row
#'   per individual, genotypes coded as 0/1/2/-9, and rownames giving individual
#'   IDs.
#' @param Par  pedigree dataframe, only genotyped parents are used.
#' @param Dups  pairs of duplicates
#'
#' @return A dataframe with columns:
#' \item{n.dam, n.sire, n.pair}{Number of dams, sires, parent-pairs successfully
#'   genotyped for the SNP}
#' \item{OHdam, OHsire}{Count of number of opposing homozygous cases}
#' \item{MEpair}{Count of Mendelian errors, includes opposing homozygous cases}
#' \item{n.dups, n.diff}{Number of duplicate pairs successfully genotyped for
#'   the SNP; number of differences}
#'
#' @seealso \code{\link{SnpStats}}.
#'
#' @keywords internal

OHperSNP <- function(GenoM, Par, Dups=NULL)
{
  GenoMx <- GenoM
  GenoMx[GenoMx<0] <- 3
  GenoMx <- rbind(GenoMx, "NA" = 3)  # add fake row for parent='NA'
  GenoMx <- GenoMx +1
  gID <- rownames(GenoM)

  # names --> rownumbers in GenoM
  PedN <- PedToNum(Par, gID, DoDummies = "no")$PedPar
  PedN[PedN==0] <- nrow(GenoMx)

  # offspring-dam-sire trio genotype combo's  ====
  Obs.OO.all <- array(0, dim=c(nrow(PedN), 3, ncol(GenoMx)))
  for (i in 1:nrow(PedN)) {
    if (PedN[i,'dam']==0 & PedN[i,'sire']==0)  next
    Obs.OO.all[i,,] <- rbind(GenoMx[i, ],
                             GenoMx[PedN[i,1], ],
                             GenoMx[PedN[i,2], ])
  }
  tbl.trio <- plyr::aaply(Obs.OO.all, 3, function(M) table(factor(M[,1], levels=1:4),
                                                          factor(M[,2], levels=1:4),
                                                          factor(M[,3], levels=1:4)))

  # count mendelian errors ===
  MER <- array(0, dim=c(4,4,4))  # offspr - mother - father
  MER[1:3,,1] <- matrix(c(0,1,2, 0,0,1, 1,0,1, 0,0,1), 3,4)  # 0/1/2/NA
  MER[1:3,,2] <- matrix(c(0,0,1, 0,0,0, 1,0,0, 0,0,0), 3,4)
  MER[1:3,,3] <- matrix(c(1,0,1, 1,0,0, 2,1,0, 1,0,0), 3,4)
  MER[1:3,,4] <- matrix(c(0,0,1, 0,0,0, 1,0,0, 0,0,0), 3,4)

  tmpdam <- apply(tbl.trio, c(1,2,3), sum)  # sum over sire genotypes (incl. missing)
  tmpsire <- apply(tbl.trio, c(1,2,4), sum)
  Counts <- data.frame(n.dam = apply(tmpdam, 1, function(M) sum(M[1:3, 1:3])),
                       n.sire = apply(tmpsire, 1, function(M) sum(M[1:3, 1:3])),
                       n.pair = apply(tbl.trio, 1, function(A) sum(A[1:3, 1:3, 1:3])),
                       OHdam = apply(tmpdam, 1, function(M) M[1,3] + M[3,1]),
                       OHsire = apply(tmpsire, 1, function(M) M[1,3] + M[3,1]),
                       MEpair = sapply(1:ncol(GenoM), function(l) sum(MER * tbl.trio[l,,,])))

  # duplicates ===
  if (!is.null(Dups)) {
    GenoNums <- setNames(seq_along(gID), gID)
    NumDups <- cbind(id1 = GenoNums[Dups[,1]],
                     id2 = GenoNums[Dups[,2]])

    OO.dup <- array(dim=c(nrow(Dups), 2, ncol(GenoMx)))
    for (i in 1:nrow(Dups)) {
      OO.dup[i,,] <- rbind(GenoMx[NumDups[i,1], ],
                           GenoMx[NumDups[i,2], ])
    }
    tbl.dup <- plyr::aaply(OO.dup, 3, function(M) table(factor(M[,1], levels=1:3),
                                                        factor(M[,2], levels=1:3)))
    Counts$n.dups <- apply(tbl.dup, 1, sum)
    Counts$n.diff <- apply(tbl.dup, 1, function(M) sum(M) - sum(diag(M)))
  }

  # out ===
  return( Counts )
}



#==============================================================================
#' @title plot SnpStats results
#'
#' @description scatter plots and histograms of allele frequency, missingness,
#'   and estimated genotyping error, across SNPs
#'
#' @param OUT  output from \code{SnpStats}
#'
#' @return plots
#'
#' @keywords internal


PlotSnpStats <- function(OUT)
{

  oldpar <- par(no.readonly = TRUE)
  oldpar <- oldpar[!names(oldpar) %in% c("pin", "fig",'plt')]
  par(mfrow=c(3,3), mai=c(.9,.8,.2,.1), xpd=NA)
  graphics::plot.new()

  OUT <- cbind(OUT,
               'MAF' = ifelse(OUT[,"AF"] <= 0.5, OUT[,"AF"], 1-OUT[,"AF"]),
               'HWE' = -log10(OUT[,"HWE.p"]))

  # top 5% worst SNPs by each metric (coloured red in scatter plots)
  q95 <- list('MAF' = OUT[,'MAF'] < stats::quantile(OUT[,'MAF'], prob=0.05, na.rm=TRUE),
              'Mis' = OUT[,"Mis"] > stats::quantile(OUT[,"Mis"], prob=0.95, na.rm=TRUE),
              'HWE' = OUT[,"HWE"] > stats::quantile(OUT[,"HWE"], prob=0.95, na.rm=TRUE))
  vars <- c('Minor allele frequency' = 'MAF',
            'Missingness' = 'Mis',
            'HWE test: -log10(p)' = 'HWE')

  for (i in seq_along(vars)) {
    for (j in seq_along(vars)) {

      par(mfg = c(i,j))  # where to put next plot

      if (i == j) {  # diagonal: histogram
        tryPlot(hist,
                OUT[,vars[i]], breaks=nrow(OUT)/5, col="grey", main="",
                xlab=names(vars[i]), cex.lab=1.3, ylab="",
                # ErrMsg = "SnpStats: Plotting area too small",
                oldpar = oldpar)

      } else if (i > j) {  # below diagonal: scatter plot
        plot(OUT[,vars[j]], OUT[,vars[i]], pch=16, cex=1.2,
             xlab=names(vars[j]), ylab=names(vars[i]), cex.lab=1.3)
        h <- max(setdiff(seq_along(vars), c(i,j)))   # 'the other one'
        points(OUT[q95[[h]], vars[j]], OUT[q95[[h]], vars[i]], pch=16, col="red")
        legend(par('usr')[1], par('usr')[4], yjust=0, # above
               legend = paste("5% worst", vars[h]), pch=16, col="red", text.col='red', box.col='red')

      } else {
        # above diagonal: empty, do nothing
      }
    }
  }

  par(oldpar)  # restore old par settings
}
