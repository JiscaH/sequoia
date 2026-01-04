#' @title Summarise Sequoia Output or Pedigree
#'
#' @description Number of assigned parents and grandparents and sibship sizes,
#'   split by genotyped, dummy, and 'observed'.
#'
#' @param SeqList the list returned by \code{\link{sequoia}}. Only elements
#'   'Pedigree' or 'PedigreePar' and 'AgePriors' are used. All ids in
#'   'PedigreePar', and only those, are presumed genotyped.
#' @param Pedigree  dataframe, pedigree with the first three columns being id -
#'   dam - sire. Column names are ignored, as are additional columns, except for
#'   columns OHdam, OHsire, MEpair, LLRdam, LLRsire, LLRpair (plotting only).
#' @param DumPrefix character vector with prefixes for dummy dams (mothers) and
#'   sires (fathers), used to distinguish between dummies and non-dummies. Will
#'   be read from \code{SeqList}'s 'Specs' if provided.
#' @param SNPd character vector with ids of SNP genotyped individuals. Only used
#'   when \code{Pedigree} is provided instead of \code{SeqList}, to distinguish
#'   between genetically assigned parents and 'observed' parents (e.g. observed
#'   in the field, or assigned previously using microsatellites). If \code{NULL}
#'   (the default), all parents are presumed observed.
#' @param Plot  show barplots and histograms of the results, as well as of the
#'   parental LLRs, Mendelian errors, and agepriors, if present.
#' @param Panels  character vector with panel(s) to plot. Choose from 'all',
#'   'G.parents' (parents of genotyped individuals), 'D.parents' (parents of
#'   dummy individuals), 'sibships' (distribution of sibship sizes), 'LLR'
#'   (log10-likelihood ratio parent/otherwise related), 'OH' (count of opposite
#'   homozygote SNPs).
#' @return A list with the following elements:
#'   \item{PedSummary}{a 2-column matrix with basic summary statistics, similar
#'   to what used to be returned by \pkg{Pedantics}' \code{pedStatSummary} (now
#'   archived on CRAN). First column refers to the complete pedigree, second
#'   column to SNP-genotyped individuals only. Maternal siblings sharing a dummy
#'   parent are counted in the 2nd column if both sibs are genotyped, but not if
#'   one of the sibs is a dummy individual.}
#'   \item{ParentCount}{an array with the number of assigned parents,
#'   split by:
#'   \itemize{
#'     \item offspringCat: Genotyped, Dummy, or Observed* (*: only when
#'     \code{Pedigree} is provided rather than \code{SeqList}, for ids which
#'     are not listed in \code{SNPd} and do not conform to \code{DumPrefix} +
#'     number (i.e. (almost) all individuals when \code{SNPd = NULL}, the
#'     default).
#'     \item offspringSex: Female, Male, Unknown, or Herm* (*: hermaphrodite,
#'     only if any individuals occur as both dam and sire). Based only on
#'     whether an individual occurs as Dam or Sire.
#'     \item parentSex: Dam or Sire
#'     \item parentCat: Genotyped, Dummy, Observed*, or None (*: as for
#'     offspringCat)
#'   }}
#'   \item{GPCount}{an array with the number of assigned grandparents,
#'   split by:
#'   \itemize{
#'     \item offspringCat: Genotyped, Dummy, Observed*, or All
#'     \item grandparent kind: maternal grandmothers (MGM),
#'   maternal grandfathers (MGF), paternal grandmothers (PGM), paternal
#'   grandfathers (PGF)
#'     \item grandparentCat: Genotyped, Dummy, Observed*, or None
#'     }}
#'   \item{SibSize}{a list with elements 'mat' (maternal half + full siblings),
#'   'pat' (paternal half + full siblings), and 'full' (full siblings). Each
#'   is a matrix with a number of rows equal to the maximum sibship size, and 3
#'   columns, splitting by the type of parent: Genotyped, Dummy, or Observed.}
#'
#' @details  The list with results is returned \code{\link{invisible}}
#'
#' @seealso \code{\link{PlotSeqSum}} to plot the output of this function;
#'   \code{\link{PlotPropAssigned}} for just a barplot of the proportion of
#'   individuals with genotyped/dummy/observed/no assigned dam/sire.
#'
#'
#' @examples
#' SummarySeq(Ped_griffin)
#' sumry_grif <- SummarySeq(SeqOUT_griffin, Panels=c("G.parents", "OH"))
#' sumry_grif$PedSummary
#'
#' # invisible results (nothing printed) when no output object is specified:
#' SummarySeq(SeqOUT_griffin, Panels=c("LLR"))
#'
#' @export

SummarySeq <- function(SeqList = NULL,
                       Pedigree = NULL,
                       DumPrefix = c("F0", "M0"),
                       SNPd = NULL,
                       Plot = TRUE,
                       Panels = "all")
{
  if(!is.null(Pedigree)) {
    PedIN <- Pedigree
  } else if ("Pedigree" %in% names(SeqList)) {
    PedIN <- SeqList[["Pedigree"]]
  } else if ("PedigreePar"  %in% names(SeqList)) {
    PedIN <- SeqList[["PedigreePar"]]
  } else if (is.data.frame(SeqList) & any(c("id", "ID") %in% names(SeqList))) {
    # makes SummarySeq() work if only pedigree is provided as (first) argument
    PedIN <- SeqList
    SeqList <- NULL
  } else {
    stop("Please provide dataframe 'Pedigree' or SeqList with element Pedigree or PedigreePar")
  }

  # check input
  Ped <- PedPolish(PedIN, ZeroToNA=TRUE, StopIfInvalid=FALSE)
    #, addParentRows=addParentRows)  <-- conflicts with getGenerations()
  Ped$Sex <- with(Ped, ifelse(id %in% dam,
                              ifelse(id %in% sire,
                                     "Herm", "Female"),
                              ifelse(id %in% sire, "Male", "Unknown")))
  if (any(Ped$Sex == "Herm")) {
    nSex <- 4
    Ped$Sex <- factor(Ped$Sex, levels=c("Female", "Male", "Unknown", "Herm"))
  } else {
    nSex <- 3
    Ped$Sex <- factor(Ped$Sex, levels=c("Female", "Male", "Unknown"))
  }

  if ("Specs" %in% names(SeqList)) {
    DPnames <- c("DummyPrefixFemale", "DummyPrefixMale", "DummyPrefixHerm")
    DumPrefix <- unlist(SeqList$Specs[intersect(names(SeqList$Specs), DPnames)])
    DumPrefix <- sapply(DumPrefix, function(x) ifelse(nchar(x)==1, paste0(x,"0"), x))
  }
  if (length(DumPrefix)==2 & any(Ped$Sex == "Herm"))
    DumPrefix <- c(DumPrefix, "H0")
  if (is.null(SNPd) & "PedigreePar" %in% names(SeqList))
    SNPd <- SeqList$PedigreePar$id

  AllPanels <- c('G.parents', 'D.parents', 'O.parents', 'sibships', 'LLR', 'OH')
  if (!all(Panels %in% c('all', AllPanels)))  stop("Invalid value for 'Panels'")
  # else issue caught by 'Plotting area too small, or other plotting problem'


  #~~~~~~~~~~~~~~~~~~~~~~~~~
  # get individual category: genotyped, dummy, observed
  getGDO <- function(id, gID = NULL, DumPrefix = NULL)
  {
    GDO <- ifelse(is.na(id),  "None",
                  ifelse(id %in% gID,  "Genotyped",
                         "Observed"))
    if (!is.null(DumPrefix)) {
      for (p in seq_along(DumPrefix)) {
        GDO[substr(id,1,nchar(DumPrefix[p])) == DumPrefix[p]] <- "Dummy"
      }
    }
    return( factor(GDO, levels=c("Genotyped", "Dummy", "Observed", "None"), ordered=TRUE) )
  }

  for (x in c("id", "dam", "sire")) {
    Ped[, paste0("GDO.", x)] <- getGDO(Ped[, x], SNPd, DumPrefix)
  }

  itypes <- c('G','D','O')
  typenames <- c("Genotyped", "Dummy", "Observed")

  #~~~~~~~~~~~~~~~~~~~~~~~~~
  ParentCount <- array(dim = c(3,nSex,2,4),
                       dimnames = list(offspringCat = itypes,
                                       offspringSex = levels(Ped$Sex),
                                       parentSex = c("Dam", "Sire"),
                                       parentCat = levels(Ped$GDO.id)))
  ParentCount[,,"Dam",] <- with(Ped, table(GDO.id, Sex, GDO.dam))[typenames,,]
  ParentCount[,,"Sire",] <- with(Ped, table(GDO.id, Sex, GDO.sire))[typenames,,]
  dimnames(ParentCount)[[1]] <- typenames

  #~~~~~~~~~~~~~~~~~~~~~~~~~
  GPX <- matrix(c("MGM", "MGF", "PGM", "PGF"), 2,2)
  PedGP <- merge(Ped,
                 setNames(Ped[, c("id","dam","sire","GDO.dam","GDO.sire")],
                          c("dam", GPX[,1], paste0("GDO.", GPX[,1]))), all.x=TRUE)
  PedGP <- merge(PedGP,
                 setNames(Ped[, c("id","dam","sire","GDO.dam","GDO.sire")],
                          c("sire",GPX[,2], paste0("GDO.", GPX[,2]))), all.x=TRUE)
  for (x in paste0("GDO.", c(GPX))) {
    PedGP[is.na(PedGP[,x]), x] <- "None"
  }
  GPCount <- array(dim=c(4,4,4),
                   dimnames = list(c(itypes, "T"),   # type of indiv
                                   c(GPX),
                                   levels(Ped$GDO.id)))  # type of GP
  for (x in c(GPX)) {
    GPCount[itypes,x,] <- table(PedGP$GDO.id, PedGP[, paste0("GDO.",x)])[typenames,]
  }
  GPCount["T",,] <- apply(GPCount[itypes,,], 2:3, sum)

  #~~~~~~~~~~~~~~~~~~~~~~~~~
  Ped$ParentPair <- with(Ped, ifelse(is.na(dam) | is.na(sire), NA,
                                     paste0(dam, "_-_", sire)))
  MaxSibSize <- c("mat" = ifelse(any(!is.na(Ped$dam)), max(table(Ped$dam)), 0),
                  "pat" = ifelse(any(!is.na(Ped$sire)), max(table(Ped$sire)), 0),
                  "full" = ifelse(any(!is.na(Ped$ParentPair)), max(table(Ped$ParentPair)), 0))
  SibSize <- list()
  for (p in c("mat", "pat")) {
    if (MaxSibSize[p] == 0)  next
    SibSize[[p]] <- matrix(NA, MaxSibSize[p], 3,
                           dimnames=list(1:MaxSibSize[p],
                                         c("Genotyped", "Dummy", "Observed")))
    for (x in c("Genotyped", "Dummy", "Observed")) {
      if (p=="mat") {
        tbl <- with(Ped, table(dam[GDO.dam==x]))
      } else if (p=="pat") {
        tbl <- with(Ped, table(sire[GDO.sire==x]))
      }
      SibSize[[p]][, x] <- table(factor(tbl, levels=1:MaxSibSize[p]))
    }
  }
  SibSize[["full"]] <- table(factor(table(Ped$ParentPair), levels=1:MaxSibSize["full"]))

  #~~~~~~~~~~~~~~~~~~~~~~~~~
  PedG <- Ped[Ped$GDO.id == "Genotyped", ]
  # sibs with dummy/observed parent still count as sibs

  if (!any(Ped$Sex == "Herm"))  {  # quicker, but not accurate with hermaphrodites
    SibCount <- matrix(NA, 3, 2,
                       dimnames = list(c("full", "mat", "pat"),
                                       c("All", "SNPd")))
    for (x in c("full", "mat", "pat")) {
      p <- c(full = "ParentPair", mat = "dam", pat = "sire")[[x]]
      for (a in 1:2) {
        if (a==1)  tbl <- table(table(Ped[, p]))
        if (a==2)  tbl <- table(table(PedG[, p]))
        SibS <- as.numeric(rownames(tbl))
        SibCount[x, a] <- sum(tbl * SibS * (SibS -1) / 2)
      }
    }
    SibCount[2:3, ] <- sweep(SibCount[2:3, ], 2, SibCount["full", ], "-")
    rownames(SibCount) <- c("full sibs", "maternal half sib", "paternal half sibs")

  } else {
    sibcats <- c("FS", "MHS", "PHS", "XHS")
    SibCount <- matrix(NA, length(sibcats), 2,
                       dimnames = list(sibcats, c("All", "SNPd")))
    for (a in 1:2) {
      if (a==1)  RelA <- GetRelM(Ped, patmat=TRUE, GenBack=1, Return="Array")
      if (a==2)  RelA <- GetRelM(PedG, patmat=TRUE, GenBack=1, Return="Array")
      for (x in sibcats) {
        SibCount[x,a] <- sum(RelA[,,x]) /2   # each pair included 2x in matrix
      }
    }
    rownames(SibCount) <- c("full sibs", "maternal half sib", "paternal half sibs",
                            "other half sibs")
  }

  dimnames(GPCount)[[1]] <- c(typenames, "All")
  dimnames(GPCount)[[2]] <- c("maternal grandmothers",
                              "maternal grandfathers",
                              "paternal grandmothers",
                              "paternal grandfathers")

  PedG$dam[PedG$GDO.dam != "Genotyped"] <- NA
  PedG$sire[PedG$GDO.sire != "Genotyped"] <- NA

  # All - SNPd
  PedSummary <- rbind("records" = c(nrow(Ped), nrow(PedG)),
                      "maternities" = c(sum(!is.na(Ped$dam)), sum(!is.na(PedG$dam))),
                      "paternities" = c(sum(!is.na(Ped$sire)), sum(!is.na(PedG$sire))),
                      SibCount,
                      apply(GPCount[c("All", "Genotyped"), ,1:3], 2:1, sum),
                      "maximum pedigree depth" = c(max(getGenerations(Ped, StopIfInvalid=FALSE)),
                                                   max(getGenerations(PedG, StopIfInvalid=FALSE))),  # ?
                      "founders" = c(sum(is.na(Ped$dam) & is.na(Ped$sire)),
                                     sum(is.na(PedG$dam) & is.na(PedG$sire)))
  )

  #~~~~~~~~~~~~~~~~~~~~~~~~~

  SummaryOUT <- list(PedSummary = PedSummary,
                     ParentCount = ParentCount,
                     GPCount = GPCount,
                     SibSize = SibSize)

  if (Plot) {
    img <- tryCatch(
      {
        suppressWarnings(PlotSeqSum(SummaryOUT, PedIN, as.character(Panels)))
      },
      error = function(e) {
        message("SummarySeq: Plotting area too small, or other plotting problem")
        return(NA)
      })
    if (!is.null(SeqList) & any(Panels=="all") & is.null(img)) {
      if (interactive()) {
        inp <- readline(prompt = "Press <Enter> to continue to next plot ...")
      }
      if ("AgePriorExtra" %in% names(SeqList)) {
        PlotAgePrior(SeqList$AgePriorExtra)
      } else {
        PlotAgePrior(SeqList$AgePriors)
      }
    }
  }
  invisible( SummaryOUT )
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Plot Summary Overview of sequoia Output
#'
#' @description visualise the numbers of assigned parents, sibship sizes, and
#'   parental LLRs
#'
#' @param SeqSum  list output from \code{\link{SummarySeq}}.
#' @param Pedigree dataframe with at least id, dam and sire in columns 1-3,
#'   respectively. If columns with parental LLRs and/or Mendelian errors are
#'   present, these will be plotted as well.
#' @param Panels  character vector with panel(s) to plot. Choose from 'all',
#'   'G.parents' (parents of genotyped individuals), 'D.parents' (parents of
#'   dummies), 'O.parents' (parents of non-genotyped non-dummies), sibships',
#'   'LLR', 'OH'.
#' @param ask  ask for user key stroke before proceeding to next plot.
#'
#' @importFrom graphics par barplot hist text axis mtext
#'
#' @examples
#' sumry <- SummarySeq(SeqOUT_griffin, Plot=FALSE)
#' PlotSeqSum(sumry, SeqOUT_griffin$Pedigree, Panels='all', ask=FALSE)
#'
#' @export

PlotSeqSum <- function(SeqSum, Pedigree=NULL, Panels="all", ask=TRUE)
{

  PanelsIN <- Panels
  AllPanels <- c('G.parents', 'D.parents', 'O.parents', 'sibships', 'LLR', 'OH')
  if (any(Panels=="all")) {
    Panels <- AllPanels
  } else if (!all(Panels %in% AllPanels)) {
    stop("Invalid value for 'Panels'")
  }
  if (!interactive())  ask <- FALSE

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # plotting specs
  oldpar <- par(no.readonly = TRUE)
  oldpar <- oldpar[!names(oldpar) %in% c("pin", "fig",'plt')]   # current plot dimensions, not setable. bug?
  par(mai=c(.9, 1.8, 1.4,.3), mfrow=c(1,1))

  col.damsire <- matrix(c("darkred", "firebrick2", "pink", "lightgrey",
                          "darkblue", "dodgerblue", "lightblue","lightgrey"),
                        4,2,
                        dimnames=list(c("G","D","O","X"), c("Dam", "Sire")))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # barplots proportion of parents genotyped / dummy / none

  barTitle <- c(G = 'No. parents and grandparents \nassigned to genotyped individuals',
                D = "No. parents assigned to \ndummy individuals",
                O = "No. parents assigned to \n non-genotyped non-dummy individuals")

  typenames <- c(G = "Genotyped", D = "Dummy", O = "Observed")

  ParCounts <- list(G = apply(SeqSum$ParentCount['Genotyped',,,], 2:3, sum),  # sum over offspringSex
                    D = SeqSum$ParentCount['Dummy',,,],   # female & male dummies separate
                    O = apply(SeqSum$ParentCount['Observed',,,], 2:3, sum))

  Nx <- sapply(ParCounts, sum)/2  # each indiv has entry for both dam & sire (incl. 'None')
  if (Nx['G'] + Nx['D'] == 0)  barTitle[['O']] <- "No. parents assigned to \n 'observed' individuals"
  for (z in c('G', 'D', 'O')) {
    if (!paste0(z,'.parents') %in% Panels)  next
    if (Nx[z] == 0) {
      Panels <- setdiff(Panels, paste0(z,'.parents'))
      if (!all(PanelsIN %in% 'all')) {
        cli::cli_alert_info('No {z}.parents panel because no {typenames[z]} individuals')  }
    }
  }

  GPCount <- list(G = SeqSum$GPCount["Genotyped",,],
                  O = SeqSum$GPCount["Observed",,])

  # prompt to click for next plot: not for first panel
  IsFirstPanel <- setNames(rep(FALSE, length(AllPanels)), AllPanels)
  IsFirstPanel[ intersect(AllPanels, Panels)[1] ] <- TRUE

  # add horizontal barplot: specify params
  barplot.ah <- function(M, ...) barplot(M, horiz = TRUE, axes = FALSE, add = TRUE,
                                         las = 1, names.arg = rep('', ncol(M)), ...)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # parents of genotyped or observed individuals

  for (z in c('G', 'O')) {
    if ( !paste0(z, '.parents') %in% Panels )  next

    if (ask & z=='O' & !IsFirstPanel['O.parents']) {
      inp <- readline(prompt = "Press <Enter> to continue to next plot ...")
    }

    bp <- barplot(t(rbind(ParCounts[[z]],
                          GPCount[[z]])[6:1,]), horiz=TRUE, las=1,
                  space=c(rep(.2,4), 1,.2),
                  main = barTitle[[z]],
                  xlab = "No. individuals",
                  names.arg=c("Pat. grandfather", "Pat. grandmother",
                              "Mat. grandfather","Mat. grandmother",
                              "Sire (father)", "Dam (mother)"))
    axis(side=3, at=c(0:10)*Nx[z]/10,
         labels=paste0(seq(0,100,10),"%"), col="darkgrey", col.axis="darkgrey")
    abline(v=c(0:10)*Nx[z]/10, col="grey", lty=3, xpd=FALSE)

    barplot.ah(t(ParCounts[[z]]["Dam",,drop=FALSE]), space=7, col=col.damsire[,1])
    barplot.ah(t(ParCounts[[z]]["Sire",,drop=FALSE]), space=5.8, col=col.damsire[,2])

    for (p in 1:2) {
      for (x in 1:4) {
        if (ParCounts[[z]][p, x]>0) {
          rot <- ParCounts[[z]][p, x]/Nx[z] < 0.05
          xx <- ifelse(x==1, ParCounts[[z]][p, x]/2,
                       ParCounts[[z]][p, x]/2 + sum(ParCounts[[z]][p, 1:(x-1)]))
          text(xx, bp[7-p], colnames(ParCounts[[z]])[x], col=ifelse(x==1, 0, 1),
               srt=ifelse(rot, 45, 0), cex=ifelse(rot, 0.8, 1))
        }
      }
    }

    barplot.ah(t(GPCount[[z]][4:1,]), space=.2)
    for (g in 1:4) {
      for (x in 1:4) {
        if (GPCount[[z]][g, x]>0) {
          rot <- GPCount[[z]][g, x]/Nx[z] < 0.05
          xx <- ifelse(x==1, GPCount[[z]][g, x]/2,
                       GPCount[[z]][g, x]/2 + sum(GPCount[[z]][g, 1:(x-1)]))
          text(xx, rev(bp[1:4])[g], colnames(GPCount[[z]])[x], col=ifelse(x==1, 0, 1),
               srt=ifelse(rot, 45, 0), cex=ifelse(rot, 0.8, 1))
        }
      }
    }
  }


  #~~~~~~~~~~~~~~~~~
  # parents of dummy individuals

  PCD <-  t(rbind(ParCounts[['D']]['Male',,],
                  ParCounts[['D']]['Female',,]))

  AnyHerm <- 'Herm' %in% dimnames(SeqSum$ParentCount)[['offspringSex']]
  if (AnyHerm)  PCD <- cbind(PCD,
                             t(ParCounts[['D']]['Herm',,]))

  if ('D.parents' %in% Panels) {

    if (ask & !IsFirstPanel['D.parents']) {
      inp <- readline(prompt = "Press <Enter> to continue to next plot ...")
    }

    bp <- barplot(PCD, horiz=TRUE, las=1, xpd=FALSE,
                  space=c(.2,.2, rep(c(1,.2), ifelse(AnyHerm, 2, 1))),
                  main= barTitle[['D']],
                  xlab = "No. individuals",
                  names.arg=rep(c("Sire", "Dam"), ifelse(AnyHerm, 3, 2)),
                  ylim = c(0, ifelse(AnyHerm, 9.4, 6.2)))
    mtext(c("Dummy females","Dummy males"), side=2, at=c(4.5,1.3), line=3, las=1, cex=1.0)
    if (AnyHerm)  mtext("Dummy\nHermaphrodites", side=2, at=7.7, line=3, las=1, cex=1.0)

    if (AnyHerm) {
      bpM <- matrix(bp, 2,3,
                    dimnames = list(c("Sire", "Dam"), c("Male", "Female","Herm")))
    } else {
      bpM <- matrix(bp, 2,2,
                    dimnames = list(c("Sire", "Dam"), c("Male", "Female")))
    }
    maxN <- max(apply(ParCounts[['D']], 1:2, sum))
    for (s in c("Male", "Female", "Herm")) {
      if (s=="Herm" & !AnyHerm)  next
      for (p in c("Dam", "Sire")) {
        barplot.ah(as.matrix(ParCounts[['D']][s,p,]), space=bpM[p,s]-0.5, col=col.damsire[,p])
        for (x in 1:4) {
          if (ParCounts[['D']][s,p,x]>0) {
            rot <- ifelse(ParCounts[['D']][s,p, x]/maxN < 0.1, TRUE, FALSE)
            xx <- ifelse(x==1, ParCounts[['D']][s,p, x]/2,
                         ParCounts[['D']][s,p, x]/2 + sum(ParCounts[['D']][s,p, 1:(x-1)]))
            text(xx, bpM[p,s], colnames(ParCounts[[z]])[x],
                 col=c(0,1,1,1)[x],
                 srt=ifelse(rot, 45, 0), cex=ifelse(rot, 0.8, 1))
          }
        }
      }
    }
  }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # hist sibship sizes
  if ('sibships' %in% Panels) {
    if (ask & !IsFirstPanel['sibships']) {
      inp <- readline(prompt = "Press <Enter> to continue to next plot ...")
    }
    np <- par(mai=c(.85, .8,1,.2), mfrow=c(1,3))
    for (p in 1:2) {
      barplot(t(SeqSum$SibSize[[p]]), col=col.damsire[,p], las=1, space=0,
              xlab="Sibship size", ylab="Count", main=paste0(c("M","P")[p], "aternal sibships"))
    }
    barplot(SeqSum$SibSize[["full"]], col="forestgreen", las=1, space=0,
            xlab="Sibship size", ylab="Count", main="Full sibships")  #  (size >1)
  }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # OH
  OHc <- c("OHdam", "OHsire", "MEpair")
  if (any(c("OH_Mother", "OH_Father", "ME_Pair") %in% names(Pedigree))) {
    these <- match(c("OH_Mother", "OH_Father", "ME_Pair"), names(Pedigree))
    names(Pedigree)[!is.na(these)] <- OHc[!is.na(these)]
  }
  if (any(OHc %in% names(Pedigree)) & 'OH' %in% Panels) {
    tbl.OH <- list()
    for (i in 1:3) {
      Pedigree[which(Pedigree[,OHc[i]]==-9), OHc[i]] <- NA
      tbl.OH[[i]] <- table(Pedigree[,OHc[i]])
    }

    if (any(!is.na(unlist(Pedigree[,OHc])))) {
      if (ask & !IsFirstPanel['OH']) {
        inp <- readline(prompt = "Press <Enter> to continue to next plot ...")
      }
      np <- par(mfrow=c(1,3), mai=c(.8,.4,.4,.1))
      Mainz.OH <- c("Offspring - dam\nopposing homozygotes", "Offspring - sire\nopposing homozygotes ",
                    "Offspring - dam - sire\nMendelian errors")
      for (i in 1:3) {
        if (!OHc[i] %in% names(Pedigree))  next
        if (all(is.na(Pedigree[,OHc[i]]))) next
        barplot(tbl.OH[[i]], ylim = c(0, max(unlist(tbl.OH))),
                space=0, las=1, col="grey", main=Mainz.OH[i], xlab="Count")
      }
    } else if (length(Panels) < length(AllPanels)) {
      cli::cli_alert_info("No `OH` panel, because OH columns are all `<NA>`")
    }
  } else if ('OH' %in% Panels & length(Panels) < length(AllPanels)) {
    cli::cli_alert_info("No `OH` panel, because no OH columns")
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # LLR
  LLRc <- c("LLRdam", "LLRsire", "LLRpair")
  if (any(LLRc %in% names(Pedigree)) & 'LLR' %in% Panels) {
    for (i in 1:3) {
      Pedigree[which(Pedigree[,LLRc[i]]==-999), LLRc[i]] <- NA
      Pedigree[which(Pedigree[,LLRc[i]]==999), LLRc[i]] <- NA
    }
    if (any(!is.na(unlist(Pedigree[,LLRc])))) {
      if (ask & !IsFirstPanel['LLR']) {
        inp <- readline(prompt = "Press <Enter> to continue to next plot ...")
      }
      np <- par(mfrow=c(1,3), mai=c(.8,.4,.4,.1))
      brks.LLR <- pretty(x=unlist(Pedigree[, LLRc]), n=50)
      Mainz.LLR <- c("LLR dam / not dam", "LLR sire / not sire", "LLR parent pair",
                     "Offspring - dam\nopposing homozygotes", "Offspring - sire\nopposing homozygotes ",
                     "Offspring - dam - sire\nMendelian errors")
      for (i in 1:3) {
        if (!LLRc[i] %in% names(Pedigree))  next
        if (all(is.na(Pedigree[,LLRc[i]]))) next
        hist(Pedigree[,LLRc[i]], breaks=brks.LLR, col="grey",
             main=Mainz.LLR[i], xlab="Log10 likelihood ratio", ylab="")
      }
    } else if (length(Panels) < length(AllPanels)) {
      cli::cli_alert_info("No `LLR` panel, because LLR columns are all `<NA>`")
    }
  } else if ('LLR' %in% Panels & length(Panels) < length(AllPanels)) {
   cli::cli_alert_info("No `LLR` panel, because no LLR columns", immediate.=TRUE)
  }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~
  par(oldpar)
}
