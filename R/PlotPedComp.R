#' @title Visualise PedCompare Output
#'
#' @description square Venn diagrams with \code{\link{PedCompare}}
#'   \code{Counts}.
#'
#' @param Counts a 7x5x2 array with counts of matches and mismatches per
#'   category (genotyped vs dummy), as returned by \code{\link{PedCompare}}.
#' @param sameSize logical, make all per-category Venn diagrams the same size
#'   \code{TRUE}, or make their size proportional to the counts (\code{FALSE},
#'   the default). If \code{TRUE}, a warning is printed at the bottom.
#'
#' @seealso \code{\link{PedCompare}}
#'
#' @examples
#' PC.g <- PedCompare(Ped1 = cbind(FieldMums_griffin, sire=NA),
#'                    Ped2 = SeqOUT_griffin$Pedigree)
#' PlotPedComp(PC.g$Counts)
#'
#' @export

PlotPedComp <- function(Counts, sameSize=FALSE) {

  oldpar <- par(no.readonly = TRUE)
  oldpar <- oldpar[!names(oldpar) %in% c("pin", "fig",'plt')]
  par(mfcol=c(1,2), mai=c(.5, .7, .7, .3), omi=c(0,0,0,0))

  COL <- c(Match="#00B000", Mismatch="#E9002D", Ped1="#009ADE", Ped2="#FFAA00")

  # Totals
  PCT <- apply(Counts["TT",,], 1, sum)
  coco <- CalcCorners(PCT)
  plot(0,0, type="n", xlim = 1.1 * range(coco[, c("xleft", "xright")]),
       ylim = c(-0.5, 1.1) * max(coco[, "ytop"]),
       axes=FALSE, xlab="", ylab="")  # , frame.plot=TRUE
  VennSquares(count=PCT, BL=c(0,0), COL=COL, withText=TRUE, withLegend=TRUE)
  mtext("TT \n dam+sire", side=3, cex=1.3)


  # per category
  cats <- c("GG", "GD", "DG", "DD")
  S <- sqrt(max(Counts["TT", "Total",]))  # distance between 'panels'

  if (sameSize) {
    RawCounts <- Counts
    Counts <- sweep(Counts, c(1,3), Counts[,"Total",], "/")
    Counts[is.na(Counts)] <- 0
    COL <- c(COL, P1only = COL[["Ped1"]], P2only = COL[["Ped2"]])
    S <- 1.2
  }

  plot(0,0, type="n", xlim = c(-0.5, 2.5)*S, ylim=c(-3.5*S, S),
       axes=FALSE, xlab="", ylab="")
  for (a in 1:4) {
    for (p in c("dam", "sire")) {
      cocos <- VennSquares(Counts[cats[a],,p],
                           BL=c(c(dam=0, sire=1.5)[p], -(a-1))*S,
                           COL=COL, withText=!sameSize, withLegend=FALSE)
      if (sameSize) {
        for (i in c("P1only", "P2only", "Mismatch", "Match")) {
          text(x = cocos$txtco[i,"x"], y = cocos$txtco[i,"y"],
               labels=RawCounts[cats[a],i,p], col=COL[i], cex=1.2)
        }
      }
    }
  }
  mtext(c("dam", "sire"), side=3, at=c(0,1.5)*S, cex=1.3)
  mtext(cats, side=2, at=seq(0.25, -2.75, by=-1)*S, las=1, cex=1.3)
  if (sameSize)  mtext("Warning: square sizes\n not proportional to N", side=1, font=3)  # italic

  par(oldpar)  # restore old par settings
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Square Venn diagram
#'
#' @description Draw Venn diagram with squares, with match/mismatch in
#'   overlapping area.
#'
#' @param count a length 5 named vector: 'Total', 'Match', 'Mismatch',
#'   'P1only', and 'P2only'.
#' @param BL  a length 2 vector with coordinates of bottom-mid of Ped1 square.
#' @param COL  a length 4 character vector with colours, named 'Match',
#'   'Mismatch', 'Ped1', 'Ped2'.
#' @param withText  logical, add count to each rectangle.
#' @param withLegend  logical, add legend at the bottom of the plot.
#'
#' @seealso \code{\link{PlotPedComp}}
#'
#' @keywords internal

VennSquares <- function(count, BL=c(0,0),
                        COL, withText=TRUE, withLegend=FALSE) {

  coco <- CalcCorners(count)
  coco[, c("xleft", "xright")] <- coco[, c("xleft", "xright")] + BL[1]
  coco[, c("ybottom","ytop")] <- coco[, c("ybottom","ytop")] + BL[2]

  txtco <- matrix(NA, 4, 2,
                  dimnames = list(c("P1only", "P2only", "Mismatch", "Match"),
                                  c("x", "y")))
  txtco["P1only","x"] <- mean(c(coco["Ped1","xleft"],
              min(coco["Match", "xleft"], coco["Ped1", "xright"], na.rm=TRUE)))
  txtco["P1only", "y"] <- mean(c(coco["Ped1", "ytop"],
              max(coco["Match", "ytop"], coco["Ped1", "ybottom"], na.rm=TRUE)))
  txtco["P2only", "x"] <- mean(c(coco["Ped2","xright"], coco["Ped1", "xright"]))
  txtco["P2only", "y"] <- mean(coco["Ped2", c("ybottom", "ytop")])
  for (i in c("Mismatch", "Match")) {
    txtco[i, "x"] <- mean(coco[i, c("xleft", "xright")])
    txtco[i, "y"] <- mean(coco[i, c("ybottom","ytop")])
  }


  for (i in c("Ped1", "Ped2", "Mismatch", "Match")) {
    do.call(rect, c(as.list(coco[i, c("xleft", "ybottom", "xright", "ytop")]),
                    col=grDevices::adjustcolor(COL[[i]], alpha.f=0.4),
                    border=COL[[i]], lwd=2))
  }

  if (withText) {
    COL <- c(COL, P1only = COL[["Ped1"]], P2only = COL[["Ped2"]])
    for (i in c("P1only", "P2only", "Mismatch", "Match")) {
      text(x = txtco[i,"x"], y = txtco[i,"y"], labels=count[i],
           col=COL[i], cex=1.2)
    }
  }

  if (withLegend)
    legend("bottom", c("Match", "Mismatch", "Ped1 only", "Ped2 only"), fill=COL,
           ncol=2, xpd=NA)

  invisible(list(coco=coco, txtco=txtco))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Corner coordinates
#'
#' @description Calculate corner coordinates for each of the four rectangles in
#'  a square Venn diagram
#'
#' @param count a length 5 named vector: 'Total', 'Match', 'Mismatch',
#'   'P1only', and 'P2only'.
#'
#' @return a 4x4 matrix with columns "xleft", "xright", "ybottom", "ytop" (as
#'   used by \code{\link{rect}}) and rows "Ped1", "Ped2", "Mismatch", "Match".
#'
#' @details the bottom-left corner of the Ped1 square is (0,0); offset is done
#'   by \code{\link{VennSquares}}. The size of the Ped1 and Ped2 squares is
#'   proportional to their count, i.e. N1 = count["Total"] -
#'   count["P2only"], and the length of each size thus proportional to the
#'   \code{sqrt} of that.
#'
#'   The x-location of the Ped2 square is a function of the amount of overlap
#'   (Match + Mismatch): if 0% overlap then coco["Ped2","xleft"] =
#'   coco["Ped1", "xright"], if 100% overlap then coco["Ped2","xright"] =
#'   coco["Ped1", "xright"]; and proportional in-between these two extremes.
#'
#'   The overlap area between Ped1 and Ped2 is split into Mismatch (bottom) and
#'   Match (top).
#'
#' @seealso \code{\link{PlotPedComp}, \link{VennSquares}}
#'
#' @keywords internal

CalcCorners <- function(count)   # IN: PedCompare() $Counts
{
  params <- c(N1 = sqrt(count[["Total"]] - count[["P2only"]]),
              N2 = sqrt(count[["Total"]] - count[["P1only"]]))  # double bracket drops name

  udiff <- function(V)  unname(diff(V))

  # distance between centroids:
  D2.range <- c(min = abs(udiff((params[c("N1","N2")])/2)),
                max = mean(params[c("N1","N2")]))
#  D2.prop <- 1 - sum(count[c("Match", "Mismatch")])/sum(count[c("Match", "Mismatch", "P2only")])
  D2.prop <- min(count[c("P1only", "P2only")])/count[["Total"]]
  params <- c(params, D2 = D2.range[["min"]] + D2.prop * udiff(D2.range))

  params <- c(params, Mis =  count[["Mismatch"]] /
                                      sum(count[c("Match", "Mismatch")])
              * min(params[c("N1","N2")]))

  coco <- matrix(NA, 4,4,
                 dimnames = list(c("Ped1", "Ped2", "Mismatch", "Match"),
                                 c("xleft", "xright", "ybottom", "ytop")))

  coco["Ped1",] <- c(0 - params["N1"]/2, 0 + params["N1"]/2,
                     0, params["N1"])
  coco["Ped2", ] <- c(params["D2"] - params["N2"]/2, params["D2"] + params["N2"]/2,
                      0, params["N2"])
  coco["Mismatch", ] <- c(params["D2"]-params["N2"]/2, params["N1"]/2,
                          0, params["Mis"])
  coco["Match", ] <- c(params["D2"]-params["N2"]/2, params["N1"]/2,
                       params["Mis"], min(params[c("N1", "N2")]))

  return( coco )
}
