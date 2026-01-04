#======================================================================
# Miscelaneous functions (and a data vector)
#======================================================================


#======================================================================
# convert named vector to 2-column data.frame
as.DF <- function(V, colnames=c("name", "x")) {
  setNames(data.frame(c1 = names(V), c2=V), colnames)
}

#======================================================================
# convert 3D array to matrix
A2M <- function(A) {
  if (dim(A)[3]!=1)  stop('function A2M only intended for arrays with dim(3)=1')
  array(A[,,1], dim=dim(A)[1:2], dimnames=dimnames(A)[1:2])
}

#======================================================================
# test if can be converted to integers/numbers ----
check.integer <- function(xx) ifelse(is.na(xx), NA,
                                     grepl("^[-]{0,1}[0-9]{1,}$", xx))


#======================================================================
# Comparison ----
eqv <- function(x, V, xNA=FALSE) {
  if (length(x)==1) {
    if (!is.na(x)) {
      y <- ifelse(!is.na(V), x==V, FALSE)
    } else if (is.na(xNA)) {
      y <- is.na(V)
    } else {
      y <- rep(xNA, length(V))
    }
  } else if (length(x)==length(V)) {
    y <- ifelse(!is.na(x) & !is.na(V),
              x==V,
              ifelse(is.na(x) & is.na(V),
                     ifelse(is.na(xNA),
                          TRUE,
                          xNA),
                     FALSE))
  } else {
    stop("unequal lengths")
  }
  y
}


#======================================================================
# Convert factor to numeric ----
FacToNum <- function(x) as.numeric(as.character(x))


#======================================================================
# inflate square matrix to larger square matrix with more IDs
inflate <- function(M, IDnew, na=NA) {
  Mnew <- matrix(na, length(IDnew), length(IDnew), dimnames=list(IDnew, IDnew))
  if (is.null(rownames(M)) & nrow(M)==ncol(M))  rownames(M) <- colnames(M)
  Mnew[rownames(M), colnames(M)] <- M
  Mnew
}


#======================================================================
# function adapted from Examples from integer {base}
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  ifelse(!is.numeric(x) | !is.finite(x),
         FALSE,
         abs(x - round(x)) < tol)
}


#======================================================================
#' @title Special Merge
#'
#' @description As regular merge, but combine data from columns with the same
#'  name.
#'
#' @param df1  first dataframe (lowest priority if \code{overwrite=TRUE}).
#' @param df2  second dataframe (highest priority if \code{overwrite=TRUE}).
#' @param by  columns used for merging, required.
#' @param overwrite  If FALSE (the default), NA's in df1 are replaced by values
#'   from df2. If TRUE, all values in df1 are overwritten by values from df2,
#'   except where df2 has NA.
#' @param ...  additional arguments to merge, such as \code{all}.
#'
#' @keywords internal 
#' @noRd

MergeFill <- function(df1, df2, by, overwrite=FALSE, ...) {
  commonNames <- names(df1)[which(colnames(df1) %in% colnames(df2))]
  commonNames <- commonNames[!commonNames %in% by]
  dfmerged <- merge(df1,df2,by=by,...)
  for(i in commonNames){
    left <- paste0(i, ".x")
    right <- paste0(i, ".y")
    if (!overwrite) {
      dfmerged[is.na(dfmerged[left]),left] <- dfmerged[is.na(dfmerged[left]),right]
    } else {
      dfmerged[!is.na(dfmerged[right]),left] <- dfmerged[!is.na(dfmerged[right]),right]
    }
    dfmerged[right]<- NULL
    colnames(dfmerged)[colnames(dfmerged) == left] <- i
  }
  dfmerged
}


#======================================================================
# make a named list, i.e. namedlist(a, b, x) i.o. list(a=a, b=b, x=x) ----
namedlist <- function(...) {
  L <- list(...)
  EnvNames <- as.character(as.list( match.call())[-1L])
  if (is.null(names(L))) {  # all elements unnamed
    names(L) <- EnvNames
  } else {  # some elements unnamed
    names(L) <- ifelse(names(L)=="", EnvNames, names(L))
  }
  return( L )
}


#======================================================================
# priority of relationships (close -> distant)
# used by GetRelM() & ComparePairs()
RelRank <- c("S", "M", "P", "MP", "O", "PO", "PO?",
               "FS","FS?", "MHS", "PHS", "HS", "HS?",
               "MGM", "MGF", "PGM", "PGF", "GP", "GO","GP?",
               "FA", "FN", "FA?", '2nd', "2nd?", "HA", "HN","HA?",
               "DFC1", "FC1", "XX?", 'Q', "Q?", "U", "X")

#======================================================================
# simpleCap ----
.simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), tolower(substring(s, 2)),
          sep = "", collapse = " ")
}


#======================================================================
# if plotting area too small: throw message & continue, instead of error stop
tryPlot <- function(FUN, ...,
#                    ErrMsg = "Plotting area too small",
                    oldpar) {
  OK <- TRUE
  img <- tryCatch(
    suppressWarnings(
      do.call(FUN, list(...))
    ),
    error = function(e) {
      message(e)  #(ErrMsg)
      return(NA)
    } )
  if (!is.null(img) && all(is.na(img))) {
    OK <- FALSE
    par(oldpar)
  }
  return(OK)
}


#======================================================================
# transform vector to matrix ----
VtoM <- function(V, nr=NULL, nc=2, Ng_odd=FALSE) {
  if(Ng_odd) {
    V <- V[1:((length(V)/nc-1)*nc)]
  }
  M <- matrix(V, length(V)/nc, nc)
  if(!is.null(nr)) M <- M[1:nr, , drop=FALSE]
  M
}


XtoM <- function(V, nr=NULL, nc=2, Ng_odd=FALSE) {
  if(Ng_odd) {
    V <- V[1 : (floor(length(V)/nc)*nc)]  # Fortran doesn't round but chops
  }
  M <- matrix(V, length(V)/nc, nc)
  if(!is.null(nr)) M <- M[1:nr, , drop=FALSE]
  M
}

#======================================================================
##' Catch *and* save both errors and warnings, and in the case of
##' a warning, also keep the computed result.
##'
##' @title tryCatch both warnings (with value) and errors
##' @param expr an \R expression to evaluate
##' @return a list with 'value' and 'warning', where
##'   'value' may be an error caught.
##' @author Martin Maechler;
##' Copyright (C) 2010-2012  The R Core Team
##'
#' @keywords internal 
#' @noRd

tryCatch.W.E <- function(expr)
{
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}


# #======================================================================
# # functions not used in current version
# #======================================================================
#
# #======================================================================
# # test if can be converted to numbers
# check.numeric <- function(xx) ifelse(is.na(xx), NA,
#                                      grepl("^[-]{0,1}[0-9]{0,}.{0,1}[0-9]{1,}$", xx))
#
#
# #======================================================================
# fc <- function(x, w=2)  formatC(x, width=w, flag="0")
#
#
# #======================================================================
# # Value Matching
# "%ina%" <- function(x, y) ifelse(!is.na(x), match(x, y, nomatch = 0) > 0, NA)
#
#
# #======================================================================
# Replace <- function(V, old, new) {
#   # base function 'replace' with match only replaces first match.
#   if (length(old) != length(new))  stop("'old' and 'new' must have same length")
#   if (!all(old %in% V))  stop("all 'old' must be in V")
#   these <- lapply(seq_along(old), function(x, y=V) which(y == old[x]))
#   newr <- rep(new, sapply(these, length))
#   replace(V, unlist(these), newr)
# }
#
#
# #======================================================================
# # table, sets UseNA to 'ifany'
# Table <- function(...) table(..., useNA="ifany")
#
#
# #======================================================================
# # create a table, and ensure that the levels TRUE, FALSE and NA are always all
# tbl.logic <- function(x) table(factor(x, levels=c(TRUE, FALSE, NA)),
#                                useNA="always")
#
#
# #======================================================================
# add GP columns to pedigree ----
# GPcols <- function(Ped) {
#   IDorder <- Ped$id   # merge() ignores sort=FALSE
#   Ped <- merge(Ped, setNames(Ped[Ped$id %in% Ped$dam, c("id","dam","sire")],
#                              c("dam", "MGM", "MGF")), all.x=TRUE)
#   Ped <- merge(Ped, setNames(Ped[Ped$id %in% Ped$sire, c("id","dam","sire")],
#                              c("sire", "PGM", "PGF")), all.x=TRUE)
#   rownames(Ped) <- Ped$id
#   ColOrder <- c("id","dam","sire","MGM", "MGF","PGM", "PGF")
#   return( Ped[IDorder, c(ColOrder, setdiff(colnames(Ped), ColOrder))] )
# }
#
#
# #======================================================================
# merge huge dataframes ----
# merge.dt <- function(df1, df2, key, quiet, ...) {
#   if (requireNamespace("data.table", quietly = TRUE)) {
#     df.12 <- as.data.frame(merge(data.table::data.table(df1, key=key),
#                    data.table::data.table(df2, key=key),
#                  ...))
#   } else {
#     df.12 <- merge(df1, df2, ...)
#     if (!quiet & (nrow(df1)>5000 | nrow(df2)>5000)) {  # fairly arbitrary
#       message("installing package 'data.table' is recommended to speed up",
#               "merging huge data.frames")
#     }
#   }
#
#   return( df.12 )
# }
