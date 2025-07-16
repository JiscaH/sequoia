#============================================================================
#' @title Compare Pairwise Relationships
#'
#' @description Compare, count and identify different types of relative pairs
#'   between two pedigrees, or within one pedigree.
#'
#' @details If \code{Pairs2} is as returned by \code{\link{GetMaybeRel}}
#'   (identified by the additional column names 'LLR' and 'OH'), these
#'   relationship categories are appended with an '?' in the output, to
#'   distinguish them from those derived from \code{Ped2}.
#'
#'   When \code{Pairs2$TopRel} contains values other than the ones listed among
#'   the return values for the combination of \code{patmat} and \code{GenBack},
#'   they are prioritised in decreasing order of factor levels, or in decreasing
#'   alphabetical order, and before the default (\code{ped2} derived) levels.
#'
#'   The matrix returned by \code{\link{DyadCompare}} [Deprecated] is a subset
#'   of the matrix returned here using default settings.
#'
#' @param Ped1 first (e.g. original/reference) pedigree, dataframe with 3
#'   columns: id-dam-sire.
#' @param Ped2 optional second (e.g. inferred) pedigree.
#' @param Pairs2 optional dataframe with as first three columns: ID1-ID2-
#'   relationship, e.g. as returned by \code{\link{GetMaybeRel}}. Column names
#'   and any additional columns are ignored. May be provided in addition to, or
#'   instead of \code{Ped2}.
#' @param GenBack  number of generations back to consider; 1 returns
#'   parent-offspring and sibling relationships, 2 also returns grandparental,
#'   avuncular and first cousins. GenBack >2 is not implemented.
#' @param patmat  logical, distinguish between paternal versus maternal relative
#'   pairs?
#' @param ExcludeDummies  logical, exclude dummy IDs from output? Individuals
#'   with e.g. the same dummy father will still be counted as paternal halfsibs.
#'   No attempt is made to match dummies in one pedigree to individuals in the
#'   other pedigree; for that use \code{\link{PedCompare}}.
#' @param  DumPrefix  character vector with the prefixes identifying dummy
#'   individuals. Use 'F0' ('M0') to avoid matching to regular individuals with
#'   IDs starting with 'F' ('M'), provided \code{Ped2} has fewer than 999 dummy
#'   females (males).
#' @param Return  return a matrix with \code{Counts} or a \code{Summary} of the
#'   number of identical relationships and mismatches per relationship, or
#'   detailed results as a 2xNxN \code{Array} or as a \code{Dataframe}.
#'   \code{All} returns a list with all four.
#' @param Pairs_suffix  symbol added to the relationship abbreviations derived
#'   from \code{Pairs2}, when both \code{Ped2} and \code{Pairs2} are provided.
#'   Can be an empty string.
#'
#' @return Depending on \code{Return}, one of the following, or a list with all:
#' \item{Counts}{(the default), a matrix with counts, with the classification in
#'   \code{Ped1} on rows and that in \code{Ped2} in columns. Counts for
#'   'symmetrical' pairs ("FS", "HS", "MHS", "PHS", "FC1", "DFC1", "U","X") are
#'   divided by two.}
#'   \item{Summary}{a matrix with one row per relationship type and four columns
#'     , named as if \code{Ped1} is the true pedigree:
#'     \describe{
#'       \item{n}{total number of pairs with that relationship in \code{Ped1},
#'         and occurring in \code{Ped2}}
#'      \item{OK}{Number of pairs with same relationship in \code{Ped2} as in
#'        \code{Ped1}}
#'      \item{hi}{Number of pairs with 'higher' relationship in \code{Ped2} as
#'        in \code{Ped1} (e.g. FS instead of HS; ranking is the order given
#'        below)}
#'      \item{lo}{Number of pairs with 'lower' relationship in \code{Ped2} as in
#'       \code{Ped1}, but not unrelated in \code{Ped2}}
#'   }}
#'   \item{Array}{a 2xNxN array (if \code{Ped2} or \code{Pairs2} is specified)
#'     or a NxN matrix , where N is the total number of individuals occurring in
#'     \code{Ped1} and/or \code{Ped2}.}
#'   \item{Dataframe}{a dataframe with \eqn{N^2} rows and four columns:
#'     \describe{
#'       \item{id.A}{First individual of the pair}
#'       \item{id.B}{Second individual of the pair}
#'       \item{RC1}{the relationship category in \code{Ped1}, as a factor with all
#'         considered categories as levels, including those with 0 count}
#'       \item{RC2}{the relationship category in \code{Ped2}}
#'     }
#'     Each pair is listed twice, e.g. once as P and once as O, or twice as FS.}
#'
#' @section Relationship abbreviations and ranking:
#' By default (\code{GenBack=1, patmat=FALSE}) the following 7 relationships are
#' distinguished:
#' \itemize{
#'    \item \strong{S}: Self (not included in \code{Counts})
#'    \item \strong{MP}: Parent
#'    \item \strong{O}: Offspring (not included in \code{Counts})
#'    \item \strong{FS}: Full sibling
#'    \item \strong{HS}: Half sibling
#'    \item \strong{U}: Unrelated, or otherwise related
#'    \item \strong{X}: Either or both individuals not occurring in both
#'      pedigrees
#' }
#' In the array and dataframe, 'MP' indicates that the second (column)
#' individual is the parent of the first (row) individual, and 'O' indicates the
#' reverse.
#'
#' When \code{GenBack=1, patmat=TRUE} the categories are (S)-M-P-(O)-FS-MHS-PHS-
#' U-X.
#'
#' When \code{GenBack=2, patmat=TRUE}, the following relationships are
#' distinguished:
#' \itemize{
#'    \item \strong{S}: Self (not included in \code{Counts})
#'    \item \strong{M}: Mother
#'    \item \strong{P}: Father
#'    \item \strong{O}: Offspring (not included in \code{Counts})
#'    \item \strong{FS}: Full sibling
#'    \item \strong{MHS}: Maternal half-sibling
#'    \item \strong{PHS}: Paternal half-sibling
#'    \item \strong{MGM}: Maternal grandmother
#'    \item \strong{MGF}: Maternal grandfather
#'    \item \strong{PGM}: Paternal grandmother
#'    \item \strong{PGF}: Paternal grandfather
#'    \item \strong{GO}: Grand-offspring (not included in \code{Counts})
#'    \item \strong{FA}: Full avuncular; maternal or paternal aunt or uncle
#'    \item \strong{HA}: Half avuncular
#'    \item \strong{FN}: Full nephew/niece (not included in \code{Counts})
#'    \item \strong{HN}: Half nephew/niece (not included in \code{Counts})
#'    \item \strong{FC1}: Full first cousin
#'    \item \strong{DFC1}: Double full first cousin
#'    \item \strong{U}: Unrelated, or otherwise related
#'    \item \strong{X}: Either or both individuals not occurring in both pedigrees
#' }
#' Note that for avuncular and cousin relationships no distinction is made
#' between paternal versus maternal, as this may differ between the two
#' individuals and would generate a large number of sub-classes. When a pair is
#' related via multiple paths, the first-listed relationship is returned. To get
#' all the different paths between a pair, use \code{\link{GetRelM}} with
#' \code{Return='Array'}.
#'
#' When \code{GenBack=2, patmat=FALSE}, MGM, MGF, PGM and PGF are combined
#' into GP, with the rest of the categories analogous to the above.
#'
#' @seealso \code{\link{PedCompare}} for individual-based comparison;
#'   \code{\link{GetRelM}} for a pairwise relationships matrix of a single
#'   pedigree; \code{\link{PlotRelPairs}} for visualisation of relationships
#'   within each pedigree.
#'
#'   To estimate P(actual relationship (Ped1) | inferred relationship (Ped2)),
#'   see examples at \code{\link{EstConf}}.
#'
#' @examples
#' PairsG <- ComparePairs(Ped_griffin, SeqOUT_griffin[["Pedigree"]],
#'                        patmat = TRUE, ExcludeDummies = TRUE, Return = "All")
#' PairsG$Counts
#'
#' # pairwise correct assignment rate:
#' PairsG$Summary[,"OK"] / PairsG$Summary[,"n"]
#'
#' # check specific pair:
#' PairsG$Array[, "i190_2010_M", "i168_2009_F"]
#' # or
#' RelDF <- PairsG$Dataframe   # for brevity
#' RelDF[RelDF$id.A=="i190_2010_M" & RelDF$id.B=="i168_2009_F", ]
#'
#' # Colony-style lists of full sib dyads & half sib dyads:
#' FullSibDyads <- with(RelDF, RelDF[Ped1 == "FS" & id.A < id.B, ])
#' HalfSibDyads <- with(RelDF, RelDF[Ped1 == "HS" & id.A < id.B, ])
#' # Use 'id.A < id.B' because each pair is listed 2x
#'
#' @export

ComparePairs <- function(Ped1 = NULL,
                         Ped2 = NULL,
                         Pairs2 = NULL,
                         GenBack = 1,
                         patmat = FALSE,
                         ExcludeDummies = TRUE,
                         DumPrefix = c("F0", "M0"),
                         Return = "Counts",
                         Pairs_suffix = '?')
{
  if(is.null(Ped1)) stop("No 'Ped1' provided")
  Return <- .simpleCap(Return)  # capitalise 1st letter
  if (!Return %in% c("Array", "Counts", "Dataframe", "Summary", "All"))
    stop("'Return' must be 'Counts', 'Array', 'Dataframe', 'Summary' or 'All'")
  if (Return == "Summary" & is.null(Ped2) & is.null(Pairs2))
    stop("Return='Summary' not availble for a single pedigree")
  # if no Ped2 & no Pairs2, 'Summary' is identical to 'Counts'.

  if (!GenBack %in% 1:2)  stop("'GenBack' must be 1 or 2")
  if (!patmat %in% c(TRUE,FALSE))  stop("'patmat' must be TRUE or FALSE")

  # check & polish pedigrees
  Ped1 <- PedPolish(Ped1, ZeroToNA=TRUE, NullOK=FALSE, StopIfInvalid=FALSE, KeepAllColumns=FALSE)
  if (!is.null(Ped2)) {
    Ped2 <- PedPolish(Ped2, ZeroToNA=TRUE, StopIfInvalid=FALSE, KeepAllColumns=FALSE)
    if (!any(Ped2$id %in% Ped1$id))  stop("no common IDs in Ped1 and Ped2")
  }

  # check & polish Pairs2
  if (!is.null(Pairs2)) {
    if (!inherits(Pairs2, c("data.frame", "matrix")))  stop("Pairs2 should be a dataframe or matrix")
    if (!any(Pairs2[,1] %in% Ped1$id) & !any(Pairs2[,2] %in% Ped1$id))  stop("no common IDs in Ped1 and Ped2")
    # is 'Pairs2' output from GetMaybeRel?
    lvls_MaybeRel <- c("PO", "FS", "HS", "GP", "FA", "2nd", "HA", "Q", "U")
    MR <- all(c("ID1", "ID2", "TopRel", "LLR", "OH") %in% names(Pairs2)) &&  all(Pairs2$TopRel %in% lvls_MaybeRel)
  } else {
    MR <- FALSE
  }

  # get relationship matrices
  RCM.1 <- GetRelM(Ped1, GenBack=GenBack, patmat=patmat, Return="Matrix")

  if (!is.null(Ped2) & is.null(Pairs2)) {
    RCM.2 <- GetRelM(Pedigree = Ped2,
                     GenBack=GenBack, patmat=patmat, Return="Matrix")
  } else if (is.null(Ped2) & !is.null(Pairs2)) {
    RCM.2 <- GetRelM(Pairs = Pairs2[, 1:3],
                     GenBack=GenBack, patmat=patmat, Return="Matrix", Pairs_suffix ='')
#    RCM.2 <- apply(RCM.2, 1, function(x) ifelse(x %in% c('U', 'S'), x, paste0(x, '?')))
  } else if (!is.null(Ped2) & !is.null(Pairs2)) {
    RCM.2 <- GetRelM(Pedigree = Ped2,
                     Pairs = Pairs2[, 1:3],
                     GenBack=GenBack, patmat=patmat, Return="Matrix",
                     Pairs_suffix = Pairs_suffix)  #ifelse(MR, '?', '_'))
  } else {
    RCM.2 <- NULL
  }

  # align the two matrices into an array
  IDs <- unique(c(colnames(RCM.1),  colnames(RCM.2)))
  RCA <- array(dim=c(2, length(IDs), length(IDs)),
               dimnames=list(c("Ped1", "Ped2"), IDs, IDs))
  RCA["Ped1", colnames(RCM.1), colnames(RCM.1)] <- RCM.1
  RCA["Ped2", colnames(RCM.2), colnames(RCM.2)] <- RCM.2
  RCA[is.na(RCA)] <- "X"

  # delete dummies (doing this earlier causes trouble with 2-generations-back GetRelM)
  if (ExcludeDummies) {
    DPnc <- nchar(DumPrefix)
    Dummies <- rep(FALSE, length(IDs))
    for (x in seq_along(DumPrefix)) {
      Dummies <- ifelse(Dummies, Dummies, substr(IDs,1,DPnc[x]) == DumPrefix[x])
    }
    if (sum(Dummies)>0) {
      RCA <- RCA[, !Dummies, !Dummies]
      IDs <- IDs[!Dummies]
    }
  }

  #@@@  return options  @@@@@@@@@@@@@@@@

  # Array ----
  if (Return == "Array") {
    if (!is.null(Ped2) | !is.null(Pairs2)) {
      return( RCA )
    } else {
      return( RCA["Ped1", , , drop=TRUE])
    }
  }

  # relationship options & order, for each GenBack/patmat combo (-> factor levels)
  lvls <- list(GB1 = list(no = c("S", "MP", "O", "FS", "HS", "U", "X"),
                          yes = c("S","M", "P", "O", "FS", "MHS", "PHS", "U", "X")),
               GB2 = list(no = c("S","MP", "O", "FS", "HS", "GP", "GO", "FA",
                                 "HA", "FN", "HN", "FC1", "DFC1", "U", "X"),
                          yes = c("S","M", "P", "O", "FS", "MHS", "PHS","MGM",
                                  "MGF", "PGM", "PGF", "GO",
                                  "FA", "FN", "HA", "HN", "DFC1", "FC1","U", "X")))

  lvls.dup.1 <- c("FS", "HS", "MHS", "PHS", "FC1", "DFC1", "U","X")   # pairs included twice in matrix with same abbreviation
  lvls.dup.2 <- lvls.dup.1

  if (!is.null(Ped2)) {  # include all possible levels in output
    lvls2.ped <- lvls[[GenBack]][[ifelse(patmat, "yes", "no")]]
  } else {
    lvls2.ped <- "X"
  }
  if (!is.null(Pairs2)) {  # only include levels that are present in Pairs2
    if (!is.null(Ped2)) {
      lvls2.pairs <- paste0(lvls_MaybeRel, Pairs_suffix)
    } else if (is.factor(Pairs2$TopRel)) {
      lvls2.pairs <- levels(Pairs2$TopRel)
    } else {
      lvls2.pairs <- unique(na.exclude(Pairs2$TopRel))
    }
    lvls.dup.2 <- unique(c(lvls.dup.2, lvls2.pairs))
  } else {
    lvls2.pairs <- "X"
  }
  lvls2 <- unique(c(lvls2.ped, lvls2.pairs))

 # if (!is.null(Ped2) & !is.null(Pairs2) & !MR) {  # different suffix used
#    RelRank <- gsub('?$', '_?', RelRank)
#  }
  lvls2 <- c(intersect(RelRank, lvls2), setdiff(lvls2, RelRank))  # sort (RelRank defined in utils.R)


  # Dataframe ----
  if (Return %in% c("Dataframe", "All")) {
    DF <- data.frame(id.A = rep(IDs, times=length(IDs)),
                     id.B = rep(IDs, each=length(IDs)),
                     Ped1 = factor(RCA["Ped1",,],
                                   levels=lvls[[GenBack]][[ifelse(patmat, "yes", "no")]]),
                     Ped2 = factor(RCA["Ped2",,], levels=lvls2))
    if (Return == "Dataframe")  return( DF )
  }

  # Counts ----
  if (Return %in% c("Counts", "Summary", "All")) {
    lvls.tbl <- list(GB1 = list(no = c("MP", "FS", "HS", "U", "X"),
                                yes = c("M", "P", "FS", "MHS", "PHS", "U", "X")),
                     GB2 = list(no = c("MP", "FS", "HS", "GP", "FA", "HA", "FC1", "DFC1", "U", "X"),
                                yes = c("M", "P", "FS", "MHS", "PHS","MGM", "MGF", "PGM", "PGF",
                                        "FA", "HA","FC1", "DFC1", "U", "X")))
    tbl <- table(Ped1 = factor(RCA["Ped1",,], levels=lvls.tbl[[GenBack]][[ifelse(patmat, "yes", "no")]]),
                 Ped2 = factor(RCA["Ped2",,], levels=lvls2[lvls2 != "S"]))

    these.dup <- rownames(tbl) %in% lvls.dup.1
    those.dup <- colnames(tbl) %in% lvls.dup.2
    tbl[these.dup, those.dup] <- tbl[these.dup, those.dup]/2
    if (Return == "Counts")  return( tbl )
  }

  # Summary ----
  if (Return %in% c("Summary", "All")) {
    if (is.null(RCM.2)) {
      ARER <- NULL  # if no Ped2 & no Pairs2, 'Summary' is identical to 'Counts'.
    } else {
      tblz <- tbl[rownames(tbl) != "X", colnames(tbl) != "X", drop=FALSE]
      rowR <- as.numeric(factor(rownames(tblz), levels=RelRank))
      colR <- as.numeric(factor(colnames(tblz), levels=RelRank))
      hi <- tblz * outer(rowR, colR, function(x,y) x > y)
      eq <- tblz * outer(rowR, colR, function(x,y) x == y)
      tblzz <- tblz
      if ("U" %in% colnames(tblz))  tblzz[,"U"] <- 0   # not if Pairs2
      lo <- tblzz * outer(rowR, colR, function(x,y) x < y)
      ARER <- cbind(n = rowSums(tblz),
                    OK = rowSums(eq),
                    lo = rowSums(lo),
                    hi = rowSums(hi))
      if (Return == "Summary")  return( ARER )
    }
  }

  # All ----
  if (Return == "All") {
    return(list(Array = RCA,
                Counts = tbl,
                Dataframe = DF,
                Summary = ARER))
  }
}


#============================================================================
#============================================================================
#' @title Compare Dyads (DEPRECATED)
#'
#' @description Count the number of half and full sibling pairs correctly and
#'   incorrectly assigned. DEPRECATED - PLEASE USE \code{\link{ComparePairs}}
#'
#' @param  Ped1 original pedigree, dataframe with 3 columns: id-dam-sire.
#' @param  Ped2 second (inferred) pedigree.
#' @param  na1  the value for missing parents in Ped1.
#'
#' @return A 3x3 table with the number of pairs assigned as full siblings (FS),
#'   half siblings (HS) or unrelated (U, including otherwise related) in the two
#'   pedigrees, with the classification in Ped1 on rows and that in Ped2 in
#'   columns.
#'
#' @seealso \code{\link{ComparePairs}} which supersedes this function;
#'   \code{\link{PedCompare}}
#'
#' @examples
#' \dontrun{
#' DyadCompare(Ped1=Ped_HSg5, Ped2=SeqOUT_HSg5$Pedigree)
#' }
#' @export

DyadCompare <- function(Ped1 = NULL,
                        Ped2 = NULL,
                        na1 = c(NA, "0"))
{
  warning("This function is deprecated, please use ComparePairs()",
          immediate.=TRUE)
  if(is.null(Ped1) || nrow(Ped1)<2) stop("No 'Ped1' provided")
  if(is.null(Ped2) || nrow(Ped2)<2) stop("No 'Ped2' provided'")
  names(Ped1)[1:3] <- c("id", "dam.1", "sire.1")
  names(Ped2)[1:3] <- c("id", "dam.2", "sire.2")
  for (i in 1:3) {
    Ped1[, i] <- as.character(Ped1[, i])
    Ped1[Ped1[, i] %in% na1, i] <- NA
  }
  for (i in 1:3) Ped2[, i] <- as.character(Ped2[, i])
  if (!any(Ped2$id %in% Ped1$id))  stop("no common IDs in Ped1 and Ped2")
  Ped1 <- PedPolish(Ped1[,1:3], ZeroToNA=TRUE, NullOK = FALSE, StopIfInvalid=FALSE, KeepAllColumns=FALSE)
  Ped2 <- PedPolish(Ped2[,1:3], ZeroToNA=TRUE, NullOK = FALSE, StopIfInvalid=FALSE, KeepAllColumns=FALSE)

  # note: each pair is counted double
  RCT <- matrix(NA, 0, 3)
  for (x in 1:nrow(Ped1)) {
    RCT <- rbind(RCT, rc(x, Ped1))
  }

  RCI <- matrix(NA, 0, 3)
  for (x in 1:nrow(Ped2)) {
    RCI <- rbind(RCI, rc(x, Ped2))
  }

  RCTI <- merge(as.data.frame(RCT, stringsAsFactors=FALSE),
                as.data.frame(RCI, stringsAsFactors=FALSE),
                by=c("id1", "id2"), all=TRUE, suffixes = c(".1", ".2"))
  RCTI <- RCTI[RCTI$id1 %in% Ped1$id & RCTI$id2 %in% Ped1$id &
                 RCTI$id1 %in% Ped2$id & RCTI$id2 %in% Ped2$id, ]
  RCTI$RC.1[is.na(RCTI$RC.1)] <- "U"
  RCTI$RC.2[is.na(RCTI$RC.2)] <- "U"
  RCTI$RC.1 <- factor(RCTI$RC.1, levels=c("FS", "HS", "U"))
  RCTI$RC.2 <- factor(RCTI$RC.2, levels=c("FS", "HS", "U"))

  tbl <- with(RCTI, table(RC.1, RC.2, useNA="ifany"))/2  # pairs included double
  tbl["U", "U"] <- nrow(Ped2) * (nrow(Ped2)-1)/2 - sum(tbl)
  tbl
  #  sweep(tbl, 1, rowSums(tbl), "/")
}

#============================================================================
#============================================================================

#' @title Find siblings

#' @param x  an ID
#' @param Ped  a pedigree with columns id - dam - sire
#'
#' @return The individuals which are full or half siblings to x, as a
#'   three-column matrix with column names id1 (x), id2 (the siblings), and
#'   RC (the relatedness category, 'FS' or 'HS').
#'
#' @keywords internal

rc <- function(x, Ped) {
  names(Ped) <- c("id", "dam", "sire")
  RelCat <- with(Ped,
                 ifelse(id == id[x], "S",
                        ifelse(eqv(dam[x],dam,FALSE) & eqv(sire[x], sire,FALSE), "FS",
                               ifelse(eqv(dam[x],dam,FALSE) |  eqv(sire[x], sire,FALSE), "HS",
                                      NA))))
  out <- cbind(id1 = Ped$id[x],
               id2 = Ped$id[!is.na(RelCat)],
               RC = stats::na.exclude(RelCat))
  out <- out[out[,"RC"] != "S", ]
  out
}

#============================================================================
#============================================================================
