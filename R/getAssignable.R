#' @title Assignability of Reference Pedigree
#'
#' @description Identify which individuals are SNP genotyped (G), and which can
#'   potentially be substituted by a dummy individual ('dummifiable', D).
#'
#' @details  Non-genotyped individuals can potentially be substituted by a dummy
#'   during pedigree reconstruction by \code{\link{sequoia}} when they have at least one genotyped
#'   offspring, and either one additional offspring (genotyped or dummy) or an
#'   genotyped/dummy parent (i.e. a grandparent to the genotyped offspring).
#'
#'   Note that this is the bare minimum requirement; e.g. grandparents are often
#'   indistinguishable from full avuncular (see \code{\link{sequoia}} and
#'   vignette for details). G-G parent-offspring pairs are only assignable if
#'   there is age information, or information from the surrounding pedigree, to
#'   tell which of the two is the parent.
#'
#'  It is assumed that all individuals in \code{SNPd} have been
#'   genotyped for a sufficient number of SNPs. To identify samples with a
#'   too-low call rate, use \code{\link{CheckGeno}}. To calculate the call rate
#'   for all samples, see the examples below.
#'
#' @param Pedigree  dataframe with columns id-dam-sire. Reference pedigree.
#' @param SNPd  character vector with ids of genotyped individuals.
#' @param minSibSize  minimum requirements to be considered dummifiable is 1
#'  genotyped offspring, and
#'   \itemize{
#'      \item '1sib1GP': at least 1 grandparent (G or D) or 1 more offspring (G
#'      or D); these are potentially assignable by \code{\link{sequoia}}
#'      \item '2sib': at least 1 more offspring (i.e. 2 siblings). Old default for
#'         \code{\link{PedCompare}}.
#'  }.
#'
#' @return The \code{Pedigree} dataframe with 3 additional columns,
#'   \code{id.cat}, \code{dam.cat} and \code{sire.cat}, with coding similar to
#'   that used by \code{\link{PedCompare}}:
#' \item{G}{Genotyped}
#' \item{D}{Dummy or 'dummifiable'}
#' \item{X}{Not genotyped and not dummifiable}
#'
#' @examples
#' PedA <- getAssignCat(Ped_HSg5, rownames(SimGeno_example))
#' tail(PedA)
#' table(PedA$dam.cat, PedA$sire.cat, useNA="ifany")
#'
#' # calculate call rate
#' \dontrun{
#' CallRates <- apply(MyGenotypes, MARGIN=1,
#'                    FUN = function(x) sum(x!=-9)) / ncol(MyGenotypes)
#' hist(CallRates, breaks=50, col="grey")
#' GoodSamples <- rownames(MyGenotypes)[ CallRates > 0.8]
#' # threshold depends on total number of SNPs, genotyping errors, proportion
#' # of candidate parents that are SNPd (sibship clustering is more prone to
#' # false positives).
#' PedA <- getAssignCat(MyOldPedigree, rownames(GoodSamples))
#' }
#
#' @export

getAssignCat <- function(Pedigree, SNPd, minSibSize = "1sib1GP") {
  # check input
  if (is.null(SNPd))  stop("Must provide 'SNPd'")
  Pedigree <- PedPolish(Pedigree, ZeroToNA=TRUE, NullOK = FALSE,
                        StopIfInvalid=FALSE, KeepAllRows = TRUE)

  if (length(intersect(Pedigree$id, SNPd)) == 0)
    stop("'Pedigree' and 'SNPd' have no IDs in common")


  IsParent <- with(Pedigree, id %in% dam | id %in% sire)
  IsFounder <- with(Pedigree, is.na(dam) & is.na(sire))
  Pedigree$id.cat <- ifelse(Pedigree$id %in% SNPd, 'G',  # genotyped
                         ifelse(IsFounder & !IsParent, 'X', # not genotyped, no parents or offspring
                                NA))  # to be determined

  # Parents with at least 1 genotyped offspring:
  Parent_G <- with(Pedigree, list(dam = na.exclude(dam[id.cat %in% 'G']),
                                  sire = na.exclude(sire[id.cat %in% 'G'])))

  for (z in 1:10) {
    Parent_todo <- with(Pedigree, list(dam = setdiff(na.exclude(dam), id[id.cat %in% c('G','D')]),
                                     sire = setdiff(na.exclude(sire), id[id.cat %in% c('G','D')])))
    # every dummy must have at least 1 genotyped offspring; remove those w only dummy offspring
    for (p in 1:2) {
      Parent_todo[[p]] <- intersect(Parent_todo[[p]], Parent_G[[p]])
    }

    n_dummies_found <- c(0,0)
    for (p in 1:2) {
      # count number of genotyped+dummy offspring for each potential dummy
      Ped_GD <- Pedigree[Pedigree$id.cat %in% c('G','D'), ]
      NumOff_GD <- table(factor(Ped_GD[,p+1], levels=Parent_todo[[p]]))

      if (minSibSize == "2sib") {
         # matcheable by PedCompare
         dummifiable_pz <- Parent_todo[[p]][NumOff_GD > 1]
      } else if (minSibSize == "1sib1GP") {
        # potentially assignable by sequoia:
        # also if 1 offspring + at least 1 genotyped/dummy parent (=sibship grandparent)
        Ped_par <- Pedigree[match(Parent_todo[[p]], Pedigree$id), ]
        HasGP_GD <- Ped_par$dam %in% Ped_GD$id | Ped_par$sire %in% Ped_GD$id
        dummifiable_pz <- Parent_todo[[p]][NumOff_GD > 1 | (NumOff_GD>0 & HasGP_GD)]
      }

      Pedigree$id.cat[ Pedigree$id %in% dummifiable_pz ] <- 'D'
      n_dummies_found[p] <- length(dummifiable_pz)
    }
    if (sum(n_dummies_found)==0)  break
  }

  Pedigree$id.cat[is.na(Pedigree$id.cat)] <- 'X'

  for (x in c("dam", "sire")) {
    Pedigree[, paste0(x, ".cat")] <- ifelse(Pedigree[[x]] %in% SNPd, 'G',  # genotyped
                                        ifelse(Pedigree[[x]] %in% Pedigree$id[Pedigree$id.cat=='D'], 'D', 'X'))
  }

  return( Pedigree )
}
