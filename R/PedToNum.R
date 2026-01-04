#=============================================================================
#' @title Turn Character Pedigree into Numeric Pedigree
#'
#' @description Genotyped individuals get rownumber in genotype matrix,
#'   non-genotyped individuals either all get an arbitrary negative number
#'   (\code{DoDummies = 'new'}) or only individuals with a dummy ID get the
#'   corresponding negative number (\code{DoDummies = 'old'}). Note that the
#'   number series will overlap for dummy males and dummy females.
#'
#' @param Pedigree dataframe with id - dam - sire. It is assumed
#'   \code{\link{PedPolish}} has been called beforehand so that column names are
#'   correct and all columns are as.character.
#' @param gID  vector with IDs of SNP-genotyped individuals.
#' @param DoDummies  'new', 'old', or 'no' (ignore all non-genotyped
#'   individuals).
#' @param DumPrefix Prefix to identify dummies when \code{DoDummies = 'old'}
#'
#' @return a list with
#'   \item{PedPar}{An nInd x 2 matrix with the numeric IDs of parents of
#'     genotyped individuals}
#'   \item{DumPar}{A matrix with parents of dummies, see
#'     \code{\link{FoldSibGPs}}}
#'   \item{Renamed}{a length-2 list (dams, sires) with each element a dataframe
#'     with columns: 'name' (original character ID), 'num' (number ID, negative)
#'     for each dummified individual}
#'   \item{Nd}{a length 2 vector, no. dummies found/created for dams and sires}
#'
#' @details  If \code{DoDummies='new'}, \code{\link{getAssignCat}} is used
#'   with \code{minSibSize ="1sib1GP"}, and any existing dummy coding is ignored
#'   (F0001, F0002 may become -3, -6). If \code{DoDummies='old'}, the existing
#'   dummy coding is respected (F0001, F0002 will become -1, -2), but other
#'   non-genotyped individuals are ignored.
#'
#' @keywords internal
#' @noRd

PedToNum <- function(Pedigree = NULL,
                     gID = NULL,
                     DoDummies = "new",
                     DumPrefix = c("F0", "M0"))
{
  if (is.null(gID))  stop("PedToNum requires gID")
  if (is.null(Pedigree)) {
    return( list(PedPar = rep(0, 2*length(gID)),
                 DumPar = rep(0, 4*as.integer(length(gID)/2)),
                 Renamed =  list(dam = data.frame(), sire=data.frame()),
                 Nd = 0) )
  } else if (!all(gID %in% Pedigree$id)) {  # including Pedigree=NULL
    Pedigree <- rbind(Pedigree,
                      data.frame(id = gID[!gID %in% Pedigree$id],
                                 dam = NA,
                                 sire = NA))
  }

  if (length(DumPrefix) > 2 & DoDummies=="old" & !all(Pedigree$id %in% gID))
    cli::cli_alert_warning(">2 DumPrefixes not supported by `PedToNum()`")

  if (!DoDummies %in% c("old", "new", "no"))
    stop("'DoDummies' must be 'old', 'new', or 'no'")
  DPnc <- nchar(DumPrefix)

  # Dummy renaming tables ----
  Renamed <- list(dam = data.frame(), sire=data.frame())
  if (DoDummies == "new") {
    PedC <- getAssignCat(Pedigree[,1:3], gID, minSibSize = "1sib1GP")
    # remove duplicates & NA's with intersect()
    dummifiable <- with(PedC, list(dam = intersect(dam, id[id.cat=='D']),
                                     sire = intersect(sire, id[id.cat=='D'])))
    for (k in 1:2) {
      Renamed[[k]] <- data.frame(name = dummifiable[[k]],
                                 num = -seq_along(dummifiable[[k]]),
                                 stringsAsFactors = FALSE)
    }
  } else if (DoDummies == "old") {
    for (k in 1:2) {
      UniqueDummies <- sort(unique(Pedigree[substr(Pedigree[,k+1], 1,
                                                   DPnc[k]) == DumPrefix[k], k+1]))
      if (length(UniqueDummies)==0)  next
      Renamed[[k]] <- data.frame(name = UniqueDummies,
                                 num = -as.numeric(substr(UniqueDummies,
                                                          DPnc[k]+1, nchar(UniqueDummies))),
                                 stringsAsFactors = FALSE)
    }
  }

  Nd <- sapply(Renamed, nrow)


  # names to numbers ----
  NumPed <- matrix(0, nrow(Pedigree), 4,
                   dimnames=list(Pedigree$id, c("id", "dam", "sire", "sex")))
  # sex used by FoldSibsGPs() to tell female/male dummies apart

  # genotyped
  GenoNums <- setNames(seq_along(gID), gID)
  for (x in 1:3) {
    NumPed[, x] <- GenoNums[Pedigree[, x]]
  }

  # dummies
  if (DoDummies %in% c("old", "new")) {
    for (k in 1:2) {    # female, male
      if (Nd[k] == 0)  next
      if (Nd[k] > 9999)  stop("Too many dummies")

      for (x in 1:3) {  # pedigree column
        if (x!=1 & x!=k+1)  next
        these <- Pedigree[,x] %in% Renamed[[k]]$name
        NumPed[these, x] <- Renamed[[k]][match(Pedigree[these,x], Renamed[[k]]$name), "num"]
        if (x==1)  NumPed[these, "sex"] <- k
      }
    }
  }
  NumPed[is.na(NumPed)] <- 0

  # fold GPs & out ----
  PedPar <- NumPed[match(gID, Pedigree$id), 2:3]
  if (DoDummies %in% c("old", "new")) {
    DumPar <- FoldSibGPs(PedNum = NumPed, Ng = length(gID), Nd = Nd)
  } else {
    DumPar <- rep(0, 4*as.integer(length(gID)/2))
  }

  return( namedlist(PedPar, DumPar, Renamed, Nd) )
}



#=============================================================================
#=============================================================================
#' @title Fold IDs of Sibship Grandparents
#'
#' @description Fold IDs of sibship grandparents into a 2 x nInd/2 x 2 array, as
#'   they are stored in Fortran, and then stretch this into a vector that can be
#'   passed to Fortran and easily be transformed back into said 3D array.
#'
#' @param PedNum pedigree, ids replaced by numbers, dummies negative.
#' @param Ng  no. genotyped indivs.
#' @param Nd length 2 vector, no. female & male dummies.
#'
#' @return An integer vector, with missing values as 0.
#'
#' @keywords internal
#' @noRd

FoldSibGPs <- function(PedNum, Ng, Nd)
{
  DumParRF <- rep(0, 4*as.integer(Ng/2))
  if (any(Nd > 0)) {
    SibshipGPs <- array(0, dim=c(2,max(Nd),2),
                        dimnames=list(c("grandma", "granddad"),
                                      1:max(Nd), c("mat", "pat")))
    for (k in 1:2) {
      if (Nd[k] == 0) next
      # pedigree subset: parents of dummy dams (k=1) or dummy sires (k=2)
      PedDum.k <- PedNum[PedNum[,"id",drop=FALSE] < 0 & PedNum[,"sex"] == k, , drop=FALSE]

      SibshipGPs[,1:Nd[k],k] <- t(PedDum.k[order(-PedDum.k[,"id"]), 2:3])
      for (s in 1:Nd[k]) {
        for (g in 1:2) {
          x <- (k-1)*2*as.integer(Ng/2) + (s-1)*2 + g
          DumParRF[x] <- SibshipGPs[g,s,k]
        }
      }
    }
  }
  return( DumParRF )
}



#=============================================================================
#=============================================================================
#' @title Change Numeric Pedigree back to Character Pedigree
#'
#' @description Reverse \code{\link{PedToNum}}, 1 column at a time.
#'
#' @param x vector with numbers.
#' @param k 1=dam, 2=sire, needed to distinguish dummy females from dummy males.
#' @param gID  vector with IDs of SNP-genotyped individuals; rownames of
#'   genotype matrix in the exact order.
#' @param DumPrefix length-2 character vector to make dummy IDs; length-3 in
#'   case of hermaphrodites.
#'
#' @return A character vector with IDs.
#'
#' @keywords internal
#' @noRd

NumToID <- function(x, k=0, gID=NULL, DumPrefix = c("F", "M"))
{
  if (length(x)==0)  return()
  if (any(is.na(x) | !is.wholenumber(x)))
    stop("x must be whole numbers, something went wrong")
  xv <- x
  xv[xv < -1e6] <- xv[xv < -1e6] + 1e6    # hermaphrodite dummy clones
  Nd.k <- ifelse(all(xv >= 0), 0, abs(min(xv, na.rm=TRUE)))

  if (Nd.k > 9999)  stop("\nMore than 9999 dummies! Cannot parse output.")
  if (Nd.k > 999 && any(nchar(DumPrefix)==1))
    cli::cli_alert_warning(c("More than 999 dummies! Please use `DummyPrefix` of >1 character to avoid ambiguity,",
                             "default is `F0` & `M0`"), wrap=TRUE)

  if (length(k)==1)  k <- rep(k, length(x))
  if (!all(k[x<0] %in% c(1,2,4)))  stop("Invalid k")

  ID <- sapply(seq_along(x), function(i) {
    ifelse(x[i] > 0,
        gID[x[i]],
        ifelse(x[i] < -1e6,
            paste0(DumPrefix[3], formatC(- (x[i] + 1e6), width=4, flag=0)),
            ifelse(x[i] < 0,
                paste0(DumPrefix[k[i]], formatC(-x[i], width=4, flag=0)),
                NA)))

  })
  return( ID )
}
