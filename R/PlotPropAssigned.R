#' @title Plot proportion of individuals that has a parent assigned
#'
#' @description For any pedigree, plot the proportion of individuals that has
#' a genotyped, dummy, observed, or no dam/sire assigned.
#'
#' @param Pedigree  dataframe where the first 3 columns are id, dam, sire.
#' @param SNPd  character vector with ids of genotyped individuals
#'   (e.g. rownames of genotype matrix).
#' @param DumPrefix character vector with prefixes for dummy dams (mothers) and
#'   sires (fathers), used to distinguish between dummies and non-dummies.
#' @param ... further arguments passed to \code{\link{barplot}}
#'
#' @details This function offers a more flexible interface to some of the plots
#'   included in \code{\link{SummarySeq}}
#'
#' @return a 2x4 matrix with counts, returned invisibly.
#'
#' @examples
#' PlotPropAssigned(SeqOUT_griffin$Pedigree, SNPd = rownames(Geno_griffin))
#'
#' @export

PlotPropAssigned <- function(Pedigree = NULL,
                               DumPrefix = c("F0", "M0"),
                               SNPd = NULL,
                               ...)   # arguments passed to barplot()
{

  Ped <- sequoia::PedPolish(Pedigree, ZeroToNA=TRUE, StopIfInvalid=FALSE,
                            addParentRows=FALSE, LoopCheck=FALSE)

  ## categorize ----
  #~~~~~~~~~~~
  getGDO <- function(id, gID = NULL, DumPrefix = NULL)
  {
    isDummy <- rep(FALSE, length(id))
    if (!is.null(DumPrefix)) {
      for (p in seq_along(DumPrefix)) {
        isDummy[ substr(id,1,nchar(DumPrefix[p])) == DumPrefix[p] ] <- TRUE
      }
    }

    GDO <- ifelse(is.na(id),  "None",
                  ifelse(id %in% gID,  "Genotyped",
                         ifelse(isDummy, "Dummy",
                                "Observed")))
    return( factor(GDO, levels=c("Genotyped", "Dummy", "Observed", "None"), ordered=TRUE) )
  }
  #~~~~~~~~~~~

  for (x in c("dam", "sire")) {
    Ped[, paste0("GDO.", x)] <- getGDO(Ped[, x], SNPd, DumPrefix)
  }


  ## count ----
  ParCounts <- rbind(Dam = c(table(Ped$GDO.dam)),
                     Sire = c(table(Ped$GDO.sire)))


  ## plot ----
  col.damsire <- matrix(c("darkred", "firebrick2", "pink", "lightgrey",
                          "darkblue", "dodgerblue", "lightblue","lightgrey"),
                        4,2,
                        dimnames=list(c("G","D","O","X"), c("Dam", "Sire")))

  Nx <- nrow(Pedigree)
  XLIM <- c(0, max(pretty(1:Nx)))
  bp <- barplot(t(ParCounts[c('Sire','Dam'),]), las=1,
                horiz = TRUE,
                xlab = "No. individuals",
                xlim = XLIM,
                #   names.arg=c("Sire (father)", "Dam (mother)"))  take names from x
                ...)

  axis(side=3, at=c(0:10)*Nx/10,
       labels=paste0(seq(0,100,10),"%"), col="darkgrey", col.axis="darkgrey")
  abline(v=c(0:10)*Nx/10, col="grey", lty=3, xpd=FALSE)

  barplot.ah <- function(M, ...) barplot(M, horiz = TRUE, axes = FALSE, add = TRUE,
                                         las = 1, names.arg = rep('', ncol(M)), ...)

  barplot.ah(t(ParCounts["Dam",,drop=FALSE]), space=1.4, col=col.damsire[,'Dam'])
  barplot.ah(t(ParCounts["Sire",,drop=FALSE]), space=0.2, col=col.damsire[,'Sire'])

  for (p in 1:2) {
    for (x in 1:4) {
      if (ParCounts[p, x]>0) {
        rot <- ParCounts[p, x]/Nx < 0.05
        xx <- ifelse(x==1, ParCounts[p, x]/2,
                     ParCounts[p, x]/2 + sum(ParCounts[p, 1:(x-1)]))
        text(xx, bp[3-p], colnames(ParCounts)[x], col=ifelse(x==1, 0, 1),
             srt=ifelse(rot, 45, 0), cex=ifelse(rot, 0.8, 1))
      }
    }
  }

  invisible(ParCounts)
}



