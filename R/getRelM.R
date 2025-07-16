#' @title Matrix with Pairwise Relationships
#'
#' @description Generate a matrix or 3D array with all pairwise relationships
#'   from a pedigree or dataframe with pairs.
#'
#' @param Pedigree  dataframe with columns id - dam - sire.
#' @param Pairs  dataframe with columns ID1 - ID2 - Rel, e.g. as returned by
#'   \code{\link{GetMaybeRel}}. Combining \code{Pedigree} and \code{Pairs} works
#'   best if the relationships are coded as listed below.
#' @param GenBack  number of generations back to consider; 1 returns
#'   parent-offspring and sibling relationships, 2 also returns grand-parental,
#'   avuncular and first cousins.
#' @param patmat  logical, distinguish between paternal versus maternal relative
#'   pairs? For avuncular pairs, the distinction is never made.
#' @param directed logical, distinguish between e.g. ID1=offspring, ID2=mother
#'   ('M') and ID1=mother, ID2=offspring ('O')? Defaults to TRUE; if FALSE both
#'   are are scored as 'PO', as are father-offspring pairs, and all
#'   grandparent-- grand-offspring pairs are scored as 'GPO', and avuncular
#'   pairs as 'FNA' and 'HNA'. Not (currently) compatible with \code{patmat}.
#'   When \code{Return}='List', each pair is included twice (as ID1-ID2 &
#'   ID2-ID1)
#' @param Return  'Matrix', 'Array', or 'List'. 'Matrix' returns an N x N matrix
#'   with the closest relationship between each pair. 'Array' returns an N x N x
#'   R array with for each of the R considered relationships whether it exists
#'   between the pair (1) or not (0). See Details below. 'List' returns a list
#'   with for each of the R considered relationships a 2-column matrix with the
#'   IDs of the pairs having such a relationship. The size of the list (in Mb)
#'   is much smaller than for the matrix or array, and this is therefore the
#'   only format suitable for pedigrees with many thousands of individuals. If
#'   \code{Pairs} is specified, the only possible return type is 'Matrix'.
#' @param Pairs_suffix  symbol added to the relationship abbreviations derived
#'   from \code{Pairs}, when both \code{Pedigree} and \code{Pairs} are provided.
#'   Can be an empty string.
#'
#' @return  If \code{Return='Matrix'}, an N x N square matrix, with N equal to
#'   the number of rows in \code{Pedigree} (after running
#'   \code{\link{PedPolish}}) or the number of unique individuals in
#'   \code{Pairs}. If \code{Return='Array'}, an N x N x R array is returned,
#'   with R, the number of different relationships, determined by \code{GenBack}
#'   and \code{patmat}.
#'
#'  The following abbreviations are used within the returned \code{Matrix}, or
#'  as names of the 3rd dimension in the \code{Array} or of the \code{List}:
#'    \item{S}{Self}
#'    \item{M}{Mother}
#'    \item{P}{Father}
#'    \item{MP}{Mother or Father (\code{patmat=FALSE})}
#'    \item{O}{Offspring}
#'    \item{FS}{Full sibling}
#'    \item{MHS}{Maternal half-sibling}
#'    \item{PHS}{Paternal half-sibling}
#'    \item{XHS}{other half-sibling (hermaphrodites)}
#'    \item{HS}{half-sibling (\code{patmat=FALSE})}
#'    \item{MGM}{Maternal grandmother}
#'    \item{MGF}{Maternal grandfather}
#'    \item{PGM}{Paternal grandmother}
#'    \item{PGF}{Paternal grandfather}
#'    \item{GP}{Grandparent (\code{patmat=FALSE})}
#'    \item{GO}{Grand-offspring}
#'    \item{FA}{Full avuncular; maternal or paternal aunt or uncle.}
#'    \item{FN}{Full nephew/niece}
#'    \item{HA}{Half avuncular}
#'    \item{HN}{Half nephew/niece}
#'    \item{DFC1}{Double full first cousin}
#'    \item{FC1}{Full first cousin}
#'    \item{U}{Unrelated (or otherwise related)}
#'    \item{X}{Unknown, e.g. when only \code{Pairs} is provided and does not
#'      include this pair}
#'
#' @details  Double relationships are ignored when \code{Return='Matrix'}, but
#'   not when \code{Return='Array'}. For example, when A and B are both
#'   mother-offspring and paternal siblings (A mated with her father to produce
#'   B), only the mother-offspring relationship will be indicated when
#'   \code{Return='Matrix'}.
#'
#'   Note that full siblings are the exception to this rule: in the \code{Array}
#'   they will be indicated as 'FS' only, and not as 'MHS' or 'PHS'. Similarly,
#'   full avuncular pairs are not indicated as 'HA'. Double half-avuncular
#'   relationships are indicated as both FA and HA.
#'
#'   When \code{Pairs} is provided, \code{GenBack} and \code{patmat} are
#'   ignored, and no check is performed if the abbreviations are compatible with
#'   other functions.
#'
#' @seealso \code{\link{ComparePairs}} for comparing pairwise relationships
#'   between two pedigrees; \code{\link{PlotRelPairs}}.
#'
#' @examples
#' Rel.griffin <- GetRelM(Ped_griffin, directed=FALSE)  # few categories
#' Rel.griffin <- GetRelM(Ped_griffin, patmat=TRUE, GenBack=2)  # many cat.
#' table(as.vector(Rel.griffin))
#' # turning matrix into vector first makes table() much faster
#' PlotRelPairs(Rel.griffin)
#'
#' @export

GetRelM <- function(Pedigree = NULL,
                    Pairs = NULL,
                    GenBack = 1,
                    patmat = FALSE,
                    directed = TRUE,
                    Return = "Matrix",
                    Pairs_suffix = '?')
{
  if (is.null(Pedigree) & is.null(Pairs))  stop("Please provide Pedigree or Pairs")
  if (!(Return %in% c("Matrix", "Array", "List")))  stop("Return must be 'Matrix', 'Array', or 'List'")
  if (!is.null(Pairs)) {
    if (Return != "Matrix")  stop("When providing Pairs, Return must be 'Matrix'")
    if (!inherits(Pairs, c("data.frame", "matrix")))  stop("Pairs should be a dataframe or matrix")
  }
  if (!directed & patmat)  stop('directed=FALSE not compatible with patmat=TRUE')

  if (!is.null(Pedigree)) {
    Pedigree <- PedPolish(Pedigree, ZeroToNA=TRUE, NullOK=FALSE)
    RelA <- GetRelA(Pedigree, GenBack = GenBack, patmat = patmat, directed=directed,
                    List = (Return == 'List'))

    if (Return == "Matrix" ) {
      rel.i <- which(RelA == 1, arr.ind=TRUE)
      rel.i <- rel.i[!duplicated(rel.i[,1:2]), ]  # 'highest' relationship only per pair
      RelNames <- dimnames(RelA)[[3]]
      RelM.ped <- matrix("U", dim(RelA)[1], dim(RelA)[2],
                         dimnames = list(Pedigree$id, Pedigree$id))
      RelM.ped[rel.i[,1:2]] <- RelNames[rel.i[,3]]
      if (is.null(Pairs))  return( RelM.ped )
    } else {
      return( RelA )
    }
  }

  if (!is.null(Pairs)) {
    Pairs <- as.data.frame(Pairs)
    names(Pairs)[1:3] <- c("ID1", "ID2", "Rel")
    for (x in 1:3)  Pairs[,x] <- as.character(Pairs[,x])

    RelM.tmp <- plyr::daply(Pairs, .variables=c("ID1", "ID2"), .fun=function(df) df$Rel,
                            .drop_o = FALSE)  # do NOT drop dimension if ID1 or ID2 is of length 1
    RelM.tmp <- A2M(RelM.tmp)  # function in Utils.R, drop 3rd dim
    # make into symmetrical matrix to get consistency in above/below diagonal:
    IDs <- unique(c(as.character(Pairs$ID1), as.character(Pairs$ID2)))
    RelM.a <- inflate(RelM.tmp, IDs)
    RelM.b <- inflate(t(RelM.tmp), IDs)
    if(any(RelM.a != RelM.b, na.rm=TRUE)) {
      stop("One or more pairs occur 2x in 'Pairs', with different relationship")
    }
    RelM.pairs <- RelM.a  # same dimensions
    RelM.pairs[,] <- "X"  # changed 2025-03-19, was 'U'
    diag(RelM.pairs) <- "S"
    RelM.pairs[!is.na(RelM.a)] <- RelM.a[!is.na(RelM.a)]
    RelM.pairs[!is.na(RelM.b)] <- RelM.b[!is.na(RelM.b)]
    if (is.null(Pedigree))  return( RelM.pairs )
  }

  # if arriving here, !is.null(Pedigree) & !is.null(Pairs)
  # assuming that Pairs comes from GetMaybeRel() --> add '?' to distinguish sources
  RelM.pairs <- apply(RelM.pairs, 1, function(x) ifelse(x %in% c('U', 'S'), x, paste0(x, '?')))

  # align IDs
  IDs <- unique(c(colnames(RelM.ped), colnames(RelM.pairs)))
  RelM.ped.i <- inflate(RelM.ped, IDs, na='X')
  RelM.pairs.i <- inflate(RelM.pairs, IDs, na='X')

  # priority of relationships (close -> distant)
  RelRank <- c("S", "M", "P", "MP", "O", 'PO', "PO?",
               "FS","FS?", "MHS", "PHS", "HS", "HS?",
               "MGM", "MGF", "PGM", "PGF", "GP", "GO", 'GPO', "GP?",
               "FA", "FN", 'FNA', "FA?", "2nd?", "HA", "HN", 'HNA', "HA?",
               "DFC1", "FC1", "XX?", "Q?", "U", "X")
  used_rels <- unique(c(RelM.ped.i, RelM.pairs.i))
  rel_lvls <- c(intersect(RelRank, used_rels), setdiff(used_rels, RelRank))

  RelM.ped.i <- factor(RelM.ped.i, levels = rel_lvls)
  RelM.pairs.i <- factor(RelM.pairs.i, levels = rel_lvls)
  RelM <- matrix(factor(pmin(as.numeric(RelM.ped.i), as.numeric(RelM.pairs.i)),
                        levels = seq_along(rel_lvls), labels=rel_lvls),
                 nrow=length(IDs), dimnames=list(IDs, IDs))

  return( RelM )
}



################################################################################
################################################################################

#' @title Array with Pairwise Relationships
#'
#' @description Generate an array indicating the relationship(s) between all
#'  pairs of individuals according to the pedigree.
#'
#' @param Ped  dataframe with columns id - dam - sire.
#' @param GenBack  number of generations back to consider; 1 returns
#'   parent-offspring and sibling relationships, 2 also returns grand-parental,
#'   avuncular and first cousins.
#' @param patmat  logical, distinguish between paternal versus maternal relative
#'   pairs? For avuncular pairs, the distinction is never made.
#' @param directed logical, distinguish between 'O' vs 'P' or group into 'PO' ?
#' @param List logical, return a list instead of the default array
#'
#' @return a 3D array indicating if the pair has the specified relationship (1)
#' or not (0). The various relationship considered are in the 3rd dimension:
#'  \item{M}{}
#'  \item{P}{}
#'  \item{FS}{full siblings, including double 'other half sibs'}
#'  \item{MS}{}
#'  \item{PS}{}
#'  \item{XS}{other sibs: mother of A is father of B, or vv}
#' etc.
#'
#' @useDynLib sequoia, .registration = TRUE
#'
#' @keywords internal

GetRelA <- function(Ped = NULL, GenBack = 1, patmat = TRUE, directed=TRUE,
                    List = FALSE)
{
  if (is.null(Ped))  stop("Please provide Pedigree")

  PedN <- PedToNum(Ped, gID = Ped[, 1], DoDummies = "no")
  nInd <- nrow(PedN$PedPar)
  nRel <- ifelse(GenBack == 1, 8, 19)

  TMP <- .Fortran(getrel,
                  nind = nInd,
                  pedrf = as.integer(PedN$PedPar),
                  nrel = as.integer(nRel),
                  relv = integer(nInd * nInd * nRel) )

  Rels <- c("S", "M", "P", "O", "FS", "MHS", "PHS", "XHS")  # XHS: 'cross'-half-sibs
  RelsB <- c("S", "MP", "O", "FS", "HS")  # patmat = FALSE
  RelsC <- c("S", "PO","FS","HS")
  if (GenBack == 2) {
    Rels <- c(Rels, "MGM", "MGF", "PGM", "PGF", "GO",
              "FA", "FN", "HA", "HN", "DFC1", "FC1")
    RelsB <- c(RelsB, "GP", "GO",
               "FA", "FN", "HA", "HN", "DFC1", "FC1")
    RelsC <- c(RelsC, "GPO", "FNA", "HNA", "DFC1", "FC1")
  }

  if (List) {  # e.g. when very large pedigree --> relA object.size too large

    IDs <- Ped[,1]
    RelL <- list()
    for (r in seq_along(Rels)) {
      tmpM <- matrix(TMP$relv[((r-1)*nInd*nInd +1) : (r*nInd*nInd)], nInd, nInd)
      tmpW <- which(tmpM == 1, arr.ind=TRUE, useNames = FALSE)
      RelL[[Rels[r]]] <- cbind(ID1 = IDs[ tmpW[,1] ],
                               ID2 = IDs[ tmpW[,2] ])
    }
    rm(TMP, tmpM, tmpW)

    if (!patmat) {   # combine maternal & paternal relatives
      RelL[['MP']] <- rbind(RelL[['M']], RelL[['P']])
      RelL[['HS']] <- rbind(RelL[['MHS']], RelL[['PHS']], RelL[['XHS']])
      if (GenBack == 2) {
        RelL[['GP']] <- rbind(RelL[['MGM']], RelL[['MGF']], RelL[['PGM']], RelL[['PGF']])
      }
      RelL <- RelL[ RelsB ]
    }
    if (!directed) {
      RelL[['PO']] <- rbind(RelL[['MP']], RelL[['O']])
      if (GenBack == 2) {
        RelL[['GPO']] <- rbind(RelL[['GP']], RelL[['GO']])
        RelL[['FNA']] <- rbind(RelL[['FA']], RelL[['FN']])
        RelL[['HNA']] <- rbind(RelL[['HA']], RelL[['HN']])
      }
      RelL <- RelL[ RelsC ]
    }

    return( RelL )

  } else {

    RelA <- array(TMP$relv, dim = c(nInd, nInd, nRel),
                  dimnames = list(ID1 = Ped[,1], ID2 = Ped[,1], Rel = Rels) )

    if (!patmat) {   # combine maternal & paternal relatives
      RelA.B <- array(0, dim = c(nInd, nInd, length(RelsB)),
                      dimnames=list(Ped$id, Ped$id, RelsB))
      RelA.B[,,"MP"] <- RelA[,,"M"] | RelA[,,"P"]
      RelA.B[,,"HS"] <- RelA[,,"MHS"] | RelA[,,"PHS"] | RelA[,,"XHS"]
      if (GenBack == 2) {
        RelA.B[,,"GP"] <- RelA[,,"MGM"] | RelA[,,"MGF"] | RelA[,,"PGM"] | RelA[,,"PGF"]
      }
      for (r in RelsB) {
        if (r %in% c("MP", "HS", "GP"))  next
        RelA.B[,,r] <- RelA[,,r]
      }
      RelA <- RelA.B
    }
    if (!directed) {
      RelA.C <- array(0, dim = c(nInd, nInd, length(RelsC)),
                      dimnames=list(Ped$id, Ped$id, RelsC))
      RelA.C[,,"PO"] <- RelA[,,"MP"] | RelA[,,"O"]
      for (r in c('S','FS','HS')) RelA.C[,,r] <- RelA[,,r]
      if (GenBack == 2) {
        RelA.C[,,'GPO'] <- RelA[,,'GP'] | RelA[,,'GO']
        RelA.C[,,'FNA'] <- RelA[,,'FA'] | RelA[,,'FN']
        RelA.C[,,'HNA'] <- RelA[,,'HA'] | RelA[,,'HN']
        for (r in c("DFC1", "FC1"))  RelA.C[,,r] <- RelA[,,r]
      }
      RelA <- RelA.C
    }

    return( RelA )
  }
}


