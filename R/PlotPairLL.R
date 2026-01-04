#' @title Plot Pair Log10-Likelihoods
#'
#' @description Colour-coded scatter plots of e.g. LLR(PO/U) against LLR(FS/U),
#'   for various relationship combinations.
#'
#' @param PairLL  dataframe, output from \code{\link{CalcPairLL}}.
#' @param combo  list with length-2 character vectors, specifying which
#'   likelihoods to plot against each other. Choose from 'PO', 'FS', 'HS', 'GP',
#'   'FA', and 'HA'. The first one gets plotted on the x-axis, the second on the
#'   y-axis. Subsequent figures will be drawn row-wise.
#' @param nrows number of rows in the figure layout. If \code{NULL}, set to
#'   \code{ceiling(length(combo)/ncols)}.
#' @param ncols number of columns in the figure layout. If both \code{nrows} and
#'   \code{ncols} are NULL, \code{ncols} is set to
#'   \code{ceiling(sqrt(length(combo)))}, and \code{nrows} will be equal to
#'   \code{ncols} or one less.
#' @param bgcol  logical, colour the upper and lower triangle background of each
#'   figure to match the specified relationship combo.
#' @param Tassign  assignment threshold, shown as grey square in bottom-left
#'   corner and a band along the diagonal.
#' @param Tfilter  filter threshold, shown as dark grey square in bottom-left.
#'
#' @details The colour of each point is determined by columns \code{focal}
#'   (outer circle) and \code{TopRel} (inner filling) of \code{PairLL}.
#'
#'   Impossible relationships (LL > 0 in \code{PairLL}) are shown as \code{-Inf}
#'   on the axes, if any are present.
#'
#' @seealso \code{\link{CalcPairLL}}.
#'
#' @importFrom graphics par layout
#'
#' @examples
#' Pairs <- data.frame(ID1 = "a01005",
#'                     ID2 = c("a00013", "a00008", "a00011", "b00001",
#'                             "b01006", "b01007", "b01013", "b01014"),
#'                     focal = rep(c("PO", "HS"), each=4))
#' PLL <- CalcPairLL(Pairs, GenoM=SimGeno_example, Plot=FALSE)
#' PlotPairLL(PLL,
#'            combo = list(c("FS", "PO"), c("HS", "FS"), c("GP", "HS"),
#'                         c("FA", "HS"), c("HA", "FA"), c("FA", "GP")),
#'            nrows = 3)
#'
#' @export

PlotPairLL <- function(PairLL,
                       combo = list(c("FS", "PO"),
                                    c("HS", "FS"),
                                    c("GP", "HS"),
                                    c("FA", "HS")),
                       nrows = NULL,
                       ncols = NULL,
                       bgcol = TRUE,
                       Tassign = 0.5,
                       Tfilter = -2.0)
{
  oldpar <- par(no.readonly = TRUE)
  oldpar <- oldpar[!names(oldpar) %in% c("pin", "fig",'plt')]

  RelNames <- c("Parent-Offspring"="PO",
                "Full Sib" = "FS",
                "Half Sib" = "HS",
                "GrandParent" = "GP",
                "Full Avuncular" = "FA",
                "Half Avuncular" = "HA")

  LL <- PairLL[, RelNames]
  LL[LL > 0] <- -Inf   # impossible.
  LLRU <- as.matrix(sweep(LL, 1, PairLL[,"U"], "-"))

  # update 'TopRel' to new Tassign
  TopRel <- PairLL$TopRel
  max2 <- function(v) max(v[v < max(v)])
  for (x in 1:nrow(PairLL)) {
    if (TopRel[x] %in% RelNames) {
      if (LLRU[x,TopRel[x]] < Tassign) {
        TopRel[x] <- "??"
      } else if ((LLRU[x,TopRel[x]] - max2(LLRU[x,])) < Tassign) {
        TopRel[x] <- "??"
      }
    }
  }



  if (!is.list(combo) || any(sapply(combo, length) != 2))
    stop("'combo' must be a list with length-2 vectors")
  if (!all(unlist(combo) %in% RelNames))
    stop("'combo' may only include PO, FS, HS, GP, FA, HA")

  # set up figure layout ----
  if (is.null(nrows) & is.null(ncols))  ncols <- ceiling(sqrt(length(combo)))
  if (is.null(nrows))  nrows <- ceiling(length(combo) / ncols)
  if (is.null(ncols))  ncols <- ceiling(length(combo) / nrows)
  if (!is.null(nrows) & !is.null(ncols)) {
    if (ncols*nrows < length(combo))  stop("won't fit: nrows * ncols < length(combo)")
  }

  mat <- matrix(c(seq_along(combo), rep(0, nrows*ncols - length(combo))),
                nrows, ncols, byrow=TRUE)
  mat <- cbind(mat, length(combo)+1)
  ly <- tryCatch(
    layout(mat, widths=c(rep.int(1, ncols), 0.5)) ,
    error = function(e) {
      message("PlotPairLL: Plotting area too small")
      return(NA)
    } )
  if (is.na(ly))  return()

  # plot ----
  # use same colours as PlotRelPairs()
  RelCol <- c(PO = "purple4", FS = "green4", HS = "cyan3",
              GP = "goldenrod2", FA = "chocolate4",
              HA= "chocolate2", U = grDevices::grey(.4),
              "2nd" = "wheat3", "??" = "wheat1")
  RelNames <- c(RelNames,
                "Unrelated" = "U",
                "2nd degree" = "2nd",
                "??" = "??")


  par(mai=c(.7,.7,.2,.2))
  for (i in seq_along(combo)) {
    OK <- tryPlot(LLRplot,
                  relx=combo[[i]][1], rely=combo[[i]][2], LLRU = LLRU,
                  fcl = PairLL$focal,
                  top = TopRel,
                  RelCol=RelCol, bgcol = bgcol,
                  Tassign=Tassign, Tfilter=Tfilter,
#                  ErrMsg = "PlotPairLL: Plotting area too small",
                  oldpar = oldpar)
    if (!OK)  return()
  }

  # legend ----
  par(mai=c(0,.3,0,0))
  OK <- tryPlot(plot,
                1,1,type="n", axes=FALSE, xlab="", ylab="",
#                ErrMsg = "PlotPairLL: Plotting area too small for legend",
                oldpar = oldpar)
  if (OK) {
    graphics::legend("center", legend=names(RelNames), fill=RelCol[RelNames],
                     title = "Relationship\nCircle: focal\nFill: TopRel",
                     title.adj=0, bty = "n", xpd=NA, cex=1.5)
  }
  par(oldpar)  # restore old settings
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Scatter Plot of Pair LLRs
#'
#' @description Plot LLR(rely/U) against LLR(relx/U), for one combination of
#'   relationships, colour coded by fcl & top.
#'
#' @param relx  relationship to plot on the x-axis. One of 'PO', 'FS', 'HS',
#'   'GP', 'FA', or 'HA'.
#' @param rely  relationship to plot on the y-axis; as \code{relx}.
#' @param LLRU  matrix with log10-likelihoods, already scaled by LL(U) for each
#'   pair.
#' @param fcl focal relationship, sets outer circle colour of points.
#' @param top  most likely relationship, sets inner filling colour of points.
#' @param RelCol  named character vector with colours to use per relationship.
#' @param bgcol  do background colour TRUE/FALSE.
#' @param Tassign  assignment threshold, shown as grey square in bottom-left
#'   corner and band along the diagonal.
#' @param Tfilter  filter threshold, shown as dark grey square in bottom-left.
#'
#' @details The background of the plot is coloured to match \code{relx} (bottom
#'   triangle) and \code{rely} (upper triangle).
#'
#' @importFrom graphics mtext abline rect polygon axis
#'
#' @keywords internal 
#' @noRd

LLRplot <- function(relx, rely, LLRU,
                    fcl, top, RelCol, bgcol,
                    Tassign=0.5, Tfilter=-2)
{

  if (length(fcl)!=nrow(LLRU) | length(top)!=nrow(LLRU))
    stop("length 'fcl' and 'top' must equal nrow(LLRU)")

  # prep ----
  LLRU <- LLRU[, c(relx, rely)]   # different xlim/ylim between plots
  if (any(is.finite(LLRU))) {
    xylim <- range(LLRU, finite=TRUE)
    xylim[1] <- min(xylim[1], 0)
  } else {
    xylim <- c(Tfilter, 3)
  }
  xylim <- range(pretty(xylim))

  if (any(!is.finite(LLRU))) {
    AnyInf <- TRUE
    xylim[1] <- xylim[1] - diff(pretty(xylim))[1] /2
    LLRU[!is.finite(LLRU)] <- xylim[1]
  } else {
    AnyInf <- FALSE
  }

  # plot ----
  ac <- function(x, a=0.7)  grDevices::adjustcolor(x, alpha.f=a)

  plot(1,1, type="n", xlim = xylim, ylim=xylim, las=1, xlab="", ylab="")
  mtext(substitute(LLR(a /U), list(a = relx)), side=1, line=3)
  mtext(substitute(LLR(b /U), list(b = rely)), side=2, line=2)
  if (bgcol) {
    polygon(x= par("usr")[c(1,1,2)], y= par("usr")[c(3,4,4)],
            col=ac(RelCol[rely], a=.2), border=NA)  # upper triangle
    polygon(x= par("usr")[c(1,2,2)], y= par("usr")[c(3,4,3)],
            col=ac(RelCol[relx], a=.2), border=NA)  # lower triangle
  }
  polygon(x= par("usr")[c(1,1,2,2)],
          y= par("usr")[c(3,3,4,4)] + c(-1,1,1,-1)*Tassign,
          col= "lightgrey", border=NA)   # band along diagonal
  rect(par("usr")[1], par("usr")[3], Tassign, Tassign,
       col="lightgrey", border=NA)   # lower-left corner
  rect(par("usr")[1], par("usr")[3], Tfilter, Tfilter,
       col="darkgrey", border=NA)
  graphics::box()
  abline(a=0, b=1, lwd=2)
  if (AnyInf) {
    axis(1, at=xylim[1], labels=expression(-infinity))
    axis(2, at=xylim[1], labels=expression(-infinity), las=1)
  }

  points(LLRU[, relx], LLRU[, rely], pch=21, lwd=1, cex=1.5,
         col = RelCol[fcl], bg = ac(RelCol[top]))
}
