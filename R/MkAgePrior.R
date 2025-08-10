# this is P(age|relationship) / P(Age), to be used in P(relationship|age) =
# P(age|relationship) * P(relationship) / P(age) where P(Age) is the emperical
# age distribution in the sample.

#======================================================================
#======================================================================

#' @title Age Priors
#'
#' @description Estimate probability ratios \eqn{P(R|A) / P(R)} for age
#'   differences A and five categories of parent-offspring and sibling
#'   relationships R.
#'
#'
#' @details \eqn{\alpha_{A,R}} is the ratio between the observed
#'   counts of pairs with age difference A and relationship R (\eqn{N_{A,R}}),
#'   and the expected counts if age and relationship were independent
#'   (\eqn{N_{.,.}*p_A*p_R}).
#'
#'   During pedigree reconstruction, \eqn{\alpha_{A,R}} are multiplied by the
#'   genetic-only \eqn{P(R|G)} to obtain a probability that the
#'   pair are relatives of type R conditional on both their age difference and
#'   their genotypes.
#'
#'   The age-difference prior is used for pairs of genotyped individuals, as
#'   well as for dummy individuals. This assumes that the propensity for a pair
#'   with a given age difference to both be sampled does not depend on their
#'   relationship, so that the ratio \eqn{P(A|R) / P(A)} does not differ between
#'   sampled and unsampled pairs.
#'
#'   For further details, see the vignette.
#
#'
#' @param Pedigree dataframe with id - dam - sire in columns 1-3, and optional
#'   column with birth years. Other columns are ignored.
#' @param LifeHistData dataframe with 3 or 5 columns: id - sex (not used) -
#'   birthyear (optional columns BY.min - BY.max - YearLast not used), with
#'   unknown birth years coded as negative numbers or NA. "Birth year" may be in
#'   any arbitrary discrete time unit relevant to the species (day, month,
#'   decade), as long as parents are never born in the same time unit as their
#'   offspring. It may include individuals not in the pedigree, and not all
#'   individuals in the pedigree need to be in LifeHistData.
#' @param MinAgeParent  minimum age of a parent, a single number (min across
#'   dams and sires) or a vector of length two (dams, sires). Defaults to 1.
#'   When there is a conflict with the minimum age in the pedigree, the pedigree
#'   takes precedent.
#' @param MaxAgeParent  maximum age of a parent, a single number (max across
#'   dams and sires) or a vector of length two (dams, sires). If NULL, it will
#'   be set to latest - earliest birth year in \code{LifeHistData}, or estimated
#'   from the pedigree if one is provided. See details below.
#' @param Discrete  discrete generations? By default (NULL), discrete
#'   generations are assumed if all parent-offspring pairs have an age
#'   difference of 1, and all siblings an age difference of 0, and there are at
#'   least 20 pairs of each category (mother, father, maternal sibling, paternal
#'   sibling). Otherwise, overlapping generations are presumed. When
#'   \code{Discrete=TRUE} (explicitly or deduced), \code{Smooth} and
#'   \code{Flatten} are always automatically set to \code{FALSE}. Use
#'   \code{Discrete=FALSE} to enforce (potential for) overlapping generations.
#' @param Flatten logical. To deal with small sample sizes for some or all
#'   relationships, calculate weighed average between the observed age
#'   difference distribution among relatives and a flat (0/1) distribution. When
#'   \code{Flatten=NULL} (the default) automatically set to TRUE when there are
#'   fewer than 20 parents with known age of either sex assigned, or fewer than
#'   20 maternal or paternal siblings with known age difference. Also advisable
#'   if the sampled relative pairs with known age difference are non-typical of
#'   the pedigree as a whole.
#' @param lambdaNW  control weighing factors when \code{Flatten=TRUE}. Weights
#'   are calculated as \eqn{W(R) = 1 - exp(-lambdaNW * N(R))}, where \eqn{N(R)}
#'   is the number of pairs with relationship R for which the age difference is
#'   known. Large values (>0.2) put strong emphasis on the pedigree, small
#'   values (<0.0001) cause the pedigree to be ignored. Default results in
#'   \eqn{W=0.5} for \eqn{N=100}.
#' @param Smooth smooth the tails of and any dips in the distribution? Sets dips
#'   (<10\% of average of neighbouring ages) to the average of the neighbouring
#'   ages, sets the age after the end (oldest observed age) to LR(end)/2, and
#'   assigns a small value (0.001) to the ages before the front (youngest
#'   observed age) and after the new end. Peaks are not smoothed out, as these
#'   are less likely to cause problems than dips, and are more likely to be
#'   genuine characteristics of the species. Is set to \code{FALSE} when
#'   generations do not overlap (\code{Discrete=TRUE}).
#' @param Plot  plot a heatmap of the results?
#' @param Return  return only a matrix with the likelihood-ratio \eqn{P(A|R) /
#'   P(A)} (\code{"LR"}) or a list including also various intermediate
#'   statistics (\code{"all"}) ?
#' @param quiet suppress messages.
#'
#'
#' @return A matrix with the probability ratio of the age difference between two
#'   individuals conditional on them being a certain type of relative
#'   (\eqn{P(A|R)}) versus being a random draw from the sample (\eqn{P(A)}).
#'   Assuming conditional independence, this equals the probability ratio of
#'   being a certain type of relative conditional on the age difference, versus
#'   being a random draw.
#'
#
#'   The matrix has one row per age difference (0 - nAgeClasses) and five
#'   columns, one for each relationship type, with abbreviations:
#'   \item{M}{Mothers}
#'   \item{P}{Fathers}
#'   \item{FS}{Full siblings}
#'   \item{MS}{Maternal half-siblings}
#'   \item{PS}{Paternal half-siblings}
#'
#'   When \code{Return}='all', a list is returned with the following elements:
#'    \item{BirthYearRange}{vector length 2}
#'    \item{MaxAgeParent}{vector length 2, see details}
#'    \item{tblA.R}{matrix with the counts per age difference (rows) /
#'      relationship (columns) combination, plus a column 'X' with age
#'      differences across all pairs of individuals}
#'   \item{PA.R}{Proportions, i.e. \code{tblA.R} divided by its \code{colSums},
#'     with full-sibling correction applied if necessary (see vignette).}
#'    \item{LR.RU.A.raw}{Proportions \code{PA.R} standardised by global age
#'      difference distribution (column 'X'); \code{LR.RU.A} prior to flattening
#'      and smoothing}
#'    \item{Weights}{vector length 4, the weights used to flatten the
#'     distributions}
#'     \item{LR.RU.A}{the ageprior, flattened and/or smoothed}
#'    \item{Specs.AP}{the names of the input \code{Pedigree} and
#'    \code{LifeHistData} (or \code{NULL}), \code{lambdaNW}, and the 'effective'
#'      settings (i.e. after any automatic update) of \code{Discrete},
#'      \code{Smooth}, and \code{Flatten}.}
#'
#'
#' @section CAUTION:
#'   The small sample correction with \code{Smooth} and/or \code{Flatten}
#'   prevents errors in one dataset, but may introduce errors in another; a
#'   single solution that fits to the wide variety of life histories and
#'   datasets is impossible. Please do inspect the matrix, e.g. with
#'   \code{PlotAgePrior}, and adjust the input parameters and/or the output
#'   matrix as necessary.
#'
#'   A few outlier birth years can heavily influence the output; these may be
#'   easiest to spot with \code{Smooth=FALSE, Flatten=FALSE}.
#'
#' @section Single cohort:
#'   When all individuals in \code{LifeHistData} have the same birth year, it is
#'   assumed that \code{Discrete=TRUE} and \code{MaxAgeParent=1}. Consequently,
#'   it is assumed there are no avuncular pairs present in the sample; cousins
#'   are considered as alternative. To enforce overlapping generations, and
#'   thereby the consideration of full- and half- avuncular relationships, set
#'   \code{MaxAgeParent} to some value greater than \eqn{1}.
#'
#'  When no birth year information is given at all, a single cohort is assumed,
#'  and the same rules apply.
#'
#'
#' @section Other time units:
#'   "Birth year" may be in any arbitrary time unit relevant to the species
#'   (day, month, decade), as long as parents are always born before their
#'   putative offspring, and never in the same time unit (e.g. parent's
#'   BirthYear= 1 (or 2001) and offspring BirthYear=5 (or 2005)). Negative
#'   numbers and NA's are interpreted as unknown, and fractional numbers are not
#'   allowed.
#'
#' @section MaxAgeParent:
#'   The maximum parental age for each sex equals the maximum of:
#'  \itemize{
#'    \item the maximum age of parents in \code{Pedigree},
#'    \item the input parameter \code{MaxAgeParent},
#'    \item the maximum range of birth years in \code{LifeHistData} (including
#'      BY.min and BY.max). Only used if both of the previous are \code{NA}, or
#'      if there are fewer than 20 parents of either sex assigned.
#'    \item 1, if \code{Discrete=TRUE} or the previous three are all \code{NA}
#'  }
#'  If the age distribution of assigned parents does not capture the maximum
#'  possible age of parents, it is advised to specify \code{MaxAgeParent} for
#'  one or both sexes. Not doing so may hinder subsequent assignment of both
#'  dummy parents and grandparents. Not compatible with \code{Smooth}. If the
#'  largest age difference in the pedigree is larger than the specified
#'  \code{MaxAgeParent}, the pedigree takes precedent (i.e. the largest of the
#'  two is used).
#'
#'  @section grandparents & avuncular
#'  The agepriors for grand-parental and avuncular pairs is calculated from
#'  these by \code{\link{sequoia}}, and included in its output as
#'  `AgePriorExtra`.
#'
#'
#' @seealso \code{\link{sequoia}} and its argument \code{args.AP},
#'   \code{\link{PlotAgePrior}} for visualisation. The age vignette gives
#'   further details, mathematical justification, and some examples.
#'
#' @examples
#' # without pedigree or lifehistdata:
#' MakeAgePrior(MaxAgeParent = c(2,3))
#' MakeAgePrior(Discrete=TRUE)
#'
#' # single cohort:
#' MakeAgePrior(LifeHistData = data.frame(ID = letters[1:5], Sex=3,
#'   BirthYear=1984))
#'
#' # overlapping generations:
#' # without pedigree: MaxAgeParent = max age difference between any pair +1
#' MakeAgePrior(LifeHistData = SeqOUT_griffin$LifeHist)
#' # with pedigree:
#' MakeAgePrior(Pedigree=Ped_griffin,
#'              LifeHistData=SeqOUT_griffin$LifeHist,
#'              Smooth=FALSE, Flatten=FALSE)
#' # with small-sample correction:
#' MakeAgePrior(Pedigree=Ped_griffin,
#'              LifeHistData=SeqOUT_griffin$LifeHist,
#'              Smooth=TRUE, Flatten=TRUE)
#'
#' # Call from sequoia() via args.AP:
#' Seq_HSg5 <- sequoia(SimGeno_example, LH_HSg5, Module="par",
#'                 args.AP=list(Discrete = TRUE),  # non-overlapping generations
#'                 CalcLLR = FALSE,   # skip time-consuming calculation of LLR's
#'                 Plot = FALSE)      # no summary plots when finished
#'
#' @export

MakeAgePrior <- function(Pedigree = NULL,
                         LifeHistData = NULL,
                         MinAgeParent = NULL,
                         MaxAgeParent = NULL,
                         Discrete = NULL,
                         Flatten = NULL,
                         lambdaNW = -log(0.5)/100,
                         Smooth = TRUE,
                         Plot = TRUE,
												 Return = "LR",
												 quiet = FALSE)
{
  call.AP <- as.list(match.call()[-1L])

  # Input check ----

	if (is.null(LifeHistData) & !is.null(Pedigree)) {
	  BY.column <- tolower(colnames(Pedigree)) %in% c("by", "birthyear")
	  if (any(BY.column)) {
	    LifeHistData <- data.frame(ID = Pedigree[,1],
                                 Sex = 3,
	                               BirthYear = Pedigree[, BY.column],
                                 stringsAsFactors = FALSE)
	  }
	}

  LifeHistData <- CheckLH(LifeHistData, gID = Pedigree[,1], sorted=FALSE)
  for (x in c("BirthYear", "BY.min", "BY.max")) {
    LifeHistData[which(LifeHistData[,x] < 0), x] <- NA
  }

	if (is.null(lambdaNW) | is.na(lambdaNW))  lambdaNW <- -log(0.5)/100  # used for Flatten
	if (!Return %in% c("LR", "all"))  stop("Return must be 'LR' or 'all'")

	if (!is.null(Discrete) && is.na(Discrete))  Discrete <- NULL
	if (!(is.null(Discrete) || Discrete %in% c(TRUE, FALSE))) {
		stop("'Discrete' must be TRUE, FALSE, or NULL")
	}
	if (!is.null(Flatten) && is.na(Flatten))  Flatten <- NULL
	if (!(is.null(Flatten) || Flatten %in% c(TRUE, FALSE))) {
		stop("'Flatten' must be TRUE, FALSE, or NULL")
	}
	if (!Smooth %in% c(TRUE, FALSE))  stop("'Smooth' must be TRUE or FALSE")

	# ~~ maximum age of parents ~~
  BYrange <- suppressWarnings(range(unlist(LifeHistData[,
                        c("BirthYear", "BY.min", "BY.max")]), na.rm=TRUE))

  if (!is.null(Discrete) && Discrete) {
    MaxAgePO <- c(1,1)
  } else {
    MaxAgePO <- rep(max(1, min(abs(diff(BYrange))+1, 99)), 2)
    # +1 allows all indivs to be siblings
    # if large age difference in LifeHistData, limit to 99
    # NOT ANYMORE: if all birth years unknown, single cohort + discrete generations assumed
  }
  # check MaxAgeParent input, change MaxAgePO if specified:
  if (!is.null(MaxAgeParent) && any(!is.na(MaxAgeParent))) {
    if (length(MaxAgeParent)==1) {
      MaxAgeParent <- rep(MaxAgeParent, 2)
    } else if (length(MaxAgeParent) > 2) {
      stop("MaxAgeParent must be NULL, a single number, or length 2 vector")
    }
    if (!is.null(Discrete) && Discrete && (MaxAgeParent[1] != MaxAgeParent[2])) {
      stop("When Discrete=TRUE, MaxAgeParent must be identical for dams & sires (?)")
    }
    if (Smooth) Smooth <- FALSE  # else maxAgePO = MaxAgeParent +2
    for (p in 1:2) {
      if (!is.na(MaxAgeParent[p]) && (MaxAgeParent[p]<=0 || !is.wholenumber(MaxAgeParent[p]))) {
        stop("'MaxAgeParent' must be a positive whole number")
      }
      if (!is.na(MaxAgeParent[p])) {
        MaxAgePO[p] <- MaxAgeParent[p]
      } # else: use default from BYrange
    }
  } else {
    MaxAgeParent <- c(NA, NA)
    if (all(is.na(BYrange)) && (is.null(Discrete) || !Discrete)) {
      stop("Must provide MaxAgeParent and/or birth years in LifeHistData when Discrete=FALSE")
    }
  }
  if (max(MaxAgePO) > 100)  stop("MaxAgePO must be smaller than 100; consider a different time unit")
	names(MaxAgePO) <- c("M", "P")

  # Check MinAgeParent
  if (is.null(MinAgeParent) || any(is.na(MinAgeParent)))  MinAgeParent <- 1
  if (length(MinAgeParent)==1) {
    MinAgePO <- rep(MinAgeParent, 2)
  } else if (length(MinAgeParent)==2) {
    MinAgePO <- MinAgeParent
  } else {
    stop("MinAgeParent must be a single number, or length 2 vector (NULL or NA also OK)")
  }
  for (p in 1:2) {
    if (!is.na(MinAgePO[p]) && (MinAgePO[p]<=0 || !is.wholenumber(MinAgePO[p]))) {
      stop("'MinAgeParent' must be a positive whole number")
    }
    if (MinAgePO[p] > MaxAgePO[p])  stop("MinAgeParent must not be larger than MaxAgeParent")
  }


  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # functions to generate & return default ageprior ----
  ReturnDefault <- function(MinP=MinAgePO, MaxP=MaxAgePO, rtrn=Return, Disc=Discrete, quietR=quiet) {
    if (is.null(Disc))  Disc <- all(MaxP == 1)
		if (!quietR)  cli::cli_alert_info(paste0("Ageprior: Flat 0/1, ",
                          ifelse(Disc, "discrete","overlapping"),
                          " generations, MaxAgeParent = ", MaxP[1],",",MaxP[2]))
    LR.RU.A.default <- MkAPdefault(MinP, MaxP, Disc)
    if (Plot)  PlotAgePrior(LR.RU.A.default)
  	if (rtrn == "all") {
	    Specs.AP <- list(Pedigree = call.AP[["Pedigree"]],
                 LifeHistData = call.AP[["LifeHistData"]],
                 Discrete = Discrete,
                 Flatten = Flatten,
                 lambdaNW = lambdaNW,
                 Smooth = Smooth)
      return( list(BirthYearRange = BYrange,
                   MaxAgeParent = MaxP,
                   tblA.R = NA,
                   RelativePairs.AgeKnown = NA,
                   Weights.AgeKnown = NA,
                   LR.RU.A.unweighed = NA,
                   LR.RU.A = LR.RU.A.default,
                   Specs.AP = Specs.AP) )
    } else {
      return( LR.RU.A.default )
    }
  }

  RR <- c("M", "P", "FS", "MS", "PS")  # relatedness categories considered


	# return if no pedigree, no birthyears, or no parents ----
	Ped <- PedPolish(Pedigree, ZeroToNA=TRUE, NullOK=TRUE)
	if (is.null(Pedigree) | all(is.na(LifeHistData$BirthYear)) |
	    sum(!is.na(Ped$dam) | !is.na(Ped$sire)) ==0) {
	  return( ReturnDefault() )
	}

  # MaxT: no. rows in tables
  MaxT <- max(MaxAgePO+1, diff(BYrange))
  if ((is.null(Discrete) || !Discrete) & Smooth)  MaxT <- MaxT +2 	# space for a smooth-tail


  # counts per age & relationship ----
  Ped.LH <- merge(Ped[, c("id", "dam", "sire")], LifeHistData[,c("ID","BirthYear")],
                  by.x="id", by.y="ID",	all.x=TRUE)  # WAS: all=TRUE
  # note: some individuals in LifeHistData may be irrelevant to pedigree of SNPd indivs

  Ped.R <- Ped.LH[! (is.na(Ped.LH$dam) & is.na(Ped.LH$sire)), ]  # individuals with at least 1 parent (quicker)
  Ped.R <- merge(PedPolish(Ped.R[,1:3]), Ped.LH[, c("id", "BirthYear")])  # BY for dropped & re-added indivs

  BYV <- setNames(Ped.R$BirthYear, Ped.R$id)
  AgeDifM <- outer(BYV, BYV, "-")
  RelL <- GetRelA(Ped.R, patmat=TRUE, GenBack=1, List = TRUE)
  tblA.R <- sapply(c("M", "P", "FS", "MHS", "PHS"),
                   function(r) table(factor(AgeDifM[ RelL[[r]] ], levels = 0:MaxT)))

  # maternal siblings = maternal half-siblings + full siblings
  tblA.R <- cbind(tblA.R,
                  "MS" = tblA.R[,"MHS"] + tblA.R[,"FS"],
                  "PS" = tblA.R[,"PHS"] + tblA.R[,"FS"])
  tblA.R <- tblA.R[, RR]

  # Reference: age difference distribution across all pairs of individuals
  AgeDifM[diag(AgeDifM)] <- NA  # self
  tblA.R <- cbind(tblA.R,
                  "X" = table(factor(AgeDifM, levels=0:MaxT)))
  tblA.R[is.na(tblA.R)] <- 0

  if (any(tblA.R["0", c("M", "P")] > 0))  stop("Some parent-offspring pairs have age difference 0")
  tblA.R["0", ] <- tblA.R["0", ] / 2   # sibs & 'X' pairs with agedif 0 counted twice


  # update MaxAgePO ----
  NAK.R <- apply(tblA.R, 2, sum, na.rm=TRUE)
  ParAgeMx <- suppressWarnings(
    apply(tblA.R[, c("M", "P")], 2, function(x) max(which(x > 0)) -1) )  # 1st row = agedif 0
  ParAgeMn <- suppressWarnings(
    apply(tblA.R[, c("M", "P")], 2, function(x) min(which(x > 0)) -1) )

  MinPairs.AgeKnown <- 20   # minimum no. mother-offspring or father-offspring pairs w known age diff
  # MaxAgeParent = user input
  DoMinWarning = FALSE
  for (p in 1:2) {
    if (MaxAgePO[p] == (diff(BYrange)+1) & is.na(MaxAgeParent[p]) & NAK.R[p] >= MinPairs.AgeKnown) {
      # ignore MaxAgePO: derived from LifeHistData
      MaxAgePO[p] <- ParAgeMx[p]
    } else {
      MaxAgePO[p] <- max(ParAgeMx[p], MaxAgePO[p], na.rm=TRUE)
    }
    if (ParAgeMn[p] < MinAgePO[p]) {   # only possible if MinAgeParent specified by user
      MinAgePO[p] <- ParAgeMn[p]
      DoMinWarning = TRUE
    }
  }

  if (any(!is.na(MaxAgeParent) & MaxAgePO > MaxAgeParent) & !quiet)
    cli::cli_alert_warning(c("In the pedigree some parents are older than `MaxAgeParent` ({MaxAgeParent}),",
                             "using new pedigree-based estimate of {MaxAgePO}"), wrap=TRUE)
  if (DoMinWarning & !quiet)
    cli::cli_alert_warning(c("In the pedigree some parents are younger than `MinAgeParent` ({MinAgeParent}),",
                             "using new pedigree-based estimate of {ParAgeMn}"), wrap=TRUE)

	# check/set Discrete ----
  if (MaxAgePO[1] == MaxAgePO[2] & all(tblA.R[- (MaxAgePO[1]+1), c("M", "P")] == 0) &  # all parents have same age
  all(tblA.R[-1, c("FS", "MS", "PS")]==0)) {  # all siblings have agedif 0
    if (all(NAK.R[c("M","P", "MS", "PS")] > MinPairs.AgeKnown)) {   # sufficient evidence for discrete gen.
      if ((is.null(Discrete) || !Discrete) & Smooth) {
        MaxT <- MaxT -2  # remove space for smooth-tail again
				tblA.R <- tblA.R[1:(MaxT +1), ]
      }
      if (is.null(Discrete))  Discrete <- TRUE
    } else if (is.null(Discrete)) {
      Discrete <- FALSE   # inconclusive: assume overlapping.
    }
  } else {  # evidence for overlapping generations
    if (!is.null(Discrete) && Discrete) {
      stop("`Discrete=TRUE`, but some parents have age >1 or some siblings have age difference >0")
    }
    Discrete <- FALSE
  }
  if (Discrete) {
    Flatten <- FALSE   # ignore user-specified
    Smooth <- FALSE
    return( ReturnDefault() )
  }


  # set Flatten ----
  if (!Discrete && any(MaxAgeParent > ParAgeMx, na.rm=TRUE)) {
    if (!quiet && !is.null(Flatten) && !Flatten) {
          cli::cli_alert_info(c("All pedigree parents younger than `MaxAgeParent`, ",
         "changed to `Flatten=TRUE` to adjust ageprior to specified `MaxAgeParent`"))
      }
      Flatten <- TRUE
  }

  if (any(NAK.R[c("M","P", "MS", "PS")] < MinPairs.AgeKnown) & !Discrete) {
    if (is.null(Flatten) || Flatten) {
      Flatten <- TRUE
    } else if (!Flatten) {
      if (any(NAK.R[c("M","P", "MS", "PS")] < 5)) {
        Flatten <- TRUE   # overrule user-specified
        cli::cli_alert_info(c("Fewer than 5 mother-offspring, father-offspring, or sibling pairs with ",
          "known age difference, changing to `Flatten=TRUE`"))
      } else {
        cli::cli_alert_warning(c("Fewer than {MinPairs.AgeKnown} mother-offspring and/or ",
        "father-offspring pair with known age difference, please consider `Flatten=TRUE`}"))
      }
    }
  } else if (is.null(Flatten)) {
    Flatten <- FALSE  # more precise AP speeds up computation
  }


	# Counts to proportions ----
  PA.R <- sweep(tblA.R, 2, NAK.R, "/")

  FSuseHS <- FALSE
  if (NAK.R["FS"] / min(NAK.R[c("MS", "PS")]) < 0.5 &   # NAK.R["FS"]==0 |
                        all(NAK.R[c("MS", "PS")]>MinPairs.AgeKnown)) {
    if (!Smooth | all(NAK.R[c("MS", "PS")] > 5*MinPairs.AgeKnown)) {
      FS.tmp <- PA.R[,"MS"] * PA.R[,"PS"]
      FS.tmp <- FS.tmp/sum(FS.tmp)
    } else {
      FS.tmp <- apply(PA.R[, c("MS", "PS")], 1, mean)
    }
    if (NAK.R["FS"]>0) {
      PA.R[,"FS"] <- (PA.R[,"FS"] + FS.tmp)/2
    } else {
      PA.R[,"FS"] <- FS.tmp
    }
    FSuseHS <- TRUE
  }

  PA.R[is.nan(PA.R)] <- 1.0

  for (r in RR) {
    if (!all(PA.R[,r] %in% c(0,1))) {
      PA.R[,r] <- PA.R[,r] / sum(PA.R[,r], na.rm=TRUE)
    }
  }
  PA.R[is.na(PA.R)] <- 0


  # prob ratio Related/Unrelated given Age difference ----
  LR.RU.A.par <- PA.R[,RR]  # same dim & dimnames
  for (r in RR) {
    if (!all(PA.R[,r] %in% c(0,1))) {
      LR.RU.A.par[,r] <- PA.R[, r] / PA.R[,"X"]
    } else {
      LR.RU.A.par[,r] <- PA.R[, r]
    }
  }
  LR.RU.A.par[!is.finite(LR.RU.A.par)] <- 0   # if PA.R[,"X"]==0


  # Flatten ----
  LR.RU.A <- LR.RU.A.par
  W.R <- 1 - exp(-lambdaNW * NAK.R[RR])
  if (Flatten) {
    # Weight of tblA.R versus flat prior, as function of sample size N
    # default: lambdaNW = -log(0.5)/100 : <50% weight if N<100, and >50% if N>100
    if (FSuseHS) {
      W.R["FS"] <- 1-exp(-lambdaNW * mean(c(NAK.R["FS"], min(NAK.R[c("MS", "PS")]))))
    }

    LR.RU.A.default <- MkAPdefault(MinP=MinAgePO,
                                   MaxP = pmax(MaxAgePO, MaxAgeParent, na.rm=TRUE),  # for input > pedigree-observed
                                   Disc=Discrete)
    if (nrow(LR.RU.A.default) < nrow(LR.RU.A.par)) {
      LR.RU.A.default <- rbind(LR.RU.A.default,
                               matrix(0, nrow(LR.RU.A.par)-nrow(LR.RU.A.default), 5))
    } else if  (nrow(LR.RU.A.default) > nrow(LR.RU.A.par)) {
      LR.RU.A.default <- LR.RU.A.default[1:nrow(LR.RU.A.par), ]
    }

    for (r in RR) {
      LR.RU.A[, r] <- W.R[r] * LR.RU.A.par[, r] + (1 - W.R[r]) * LR.RU.A.default[, r]
    }
  }



  # Smooth ----
  if (Smooth) {
    LR.RU.A <- apply(LR.RU.A, 2, SmoothAP, tiny = 0.001)
    MaxAgePO <- MaxAgePO +2   # smoothed tail
  }
  LR.RU.A[1, c("M", "P")] <- 0
  LR.RU.A[is.na(LR.RU.A)] <- 0
  LR.RU.A <- round(LR.RU.A,3)


  # Out ----

  # safety check
  if (!all(apply(LR.RU.A, 2, function(x) any(x > 0))))
    stop("AgePriors error: some relationships are impossible for all age differences")

  Specs.AP <- list(Pedigree = call.AP[["Pedigree"]],
                   LifeHistData = call.AP[["LifeHistData"]],
                   Discrete = Discrete,
                   Flatten = Flatten,
                   lambdaNW = lambdaNW,
                   Smooth = Smooth)

  OUT <- list(BirthYearRange = BYrange,
              MaxAgeParent = MaxAgePO,
              tblA.R = tblA.R,
              PA.R = PA.R,
              LR.RU.A.raw = round(LR.RU.A.par,3),
              Weights = round(W.R,4),
              LR.RU.A = LR.RU.A,
              Specs.AP = Specs.AP)
  if (Plot) {
    PlotAgePrior( AP = OUT[["LR.RU.A"]] )
  }
	if (!quiet)  cli::cli_alert_info(paste0("Ageprior: Pedigree-based, ",
	                     ifelse(Discrete, "discrete ","overlapping "), "generations",
	                     ifelse(Flatten, ", flattened",""),
	                     ifelse(Smooth, ", smoothed",""),
	                     ", MaxAgeParent = ", MaxAgePO[1], ",", MaxAgePO[2]))
  utils::flush.console()
  if (Return == "all") {
    return( OUT )
  } else {
    MaxA <- min(max(which(apply(LR.RU.A, 1, function(x) any(x>0)))) +1, nrow(LR.RU.A))
    return( LR.RU.A[1:MaxA, ])
  }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Make default 1/0 ageprior
#'
#' @description Create ageprior matrix based on min and max age of parents.
#'
#' @param MinP  Minimum age of dams, sires
#' @param MaxP  Maximum age of dams, sires
#' @param Disc  Discrete generations? TRUE/FALSE/NULL
#'
#' @keywords internal

MkAPdefault <- function(MinP, MaxP, Disc) {
  RR <- c("M", "P", "FS", "MS", "PS")
  AP <- matrix(1, max(MaxP)+2, 5,
               dimnames=list(0:(max(MaxP)+1), RR))
  if(!is.null(Disc) && Disc) {  # always: MaxP[1]==MaxP[2]
    AP[,] <- 0
    AP[MaxP[1]+1, c("M","P")] <- 1
    AP[1, c("FS","MS","PS")] <- 1

  } else {
    AP[,] <- 1
    AP[1:MinP[1], "M"] <- 0
    AP[1:MinP[2], "P"] <- 0
    AP[(MaxP[1]+2):nrow(AP), "M"] <- 0
    AP[(MaxP[2]+2):nrow(AP), "P"] <- 0
    AP[(MaxP[1]-MinP[1]+2):nrow(AP), c("FS", "MS")] <- 0
    AP[(MaxP[2]-MinP[2]+2):nrow(AP), c("FS", "PS")] <- 0
  }
  return( AP )
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Smooth out dips in ageprior matrix
#'
#' @description ...
#'
#' @param V  column in ageprior matrix (vector); strictly positive
#' @param tiny  smallest non-zero value in V
#'
#' @details Sets dips (<10\% of average of neighbouring ages) to the average of
#'   the neighbouring ages, sets the age after the end (oldest observed age) to
#'   LR(end)/2, and assigns a small value (0.001) to the ages before the front
#'   (youngest observed age) and after the new end. Peaks are not smoothed out,
#'   as these are less likely to cause problems than dips, and are more likely
#'   to be genuine characteristics of the species.
#'
#' @keywords internal

SmoothAP <- function(V, tiny=0.001) {
  Front <- max(1, min(which(V > tiny)), na.rm=TRUE)
  End <- min(max(which(V > tiny)), length(V), na.rm=TRUE)
  lowMid <- rep(FALSE, length(V))
  if (End - Front >= 2) {
    for (i in c((Front+1) : (End-1))) {
      lowMid[i] <- ifelse(V[i] <= tiny, TRUE,
                      ifelse(V[i] / ((V[i-1] + V[i+1])/2) < 0.1, TRUE, FALSE))
    }
  }
  W <- V
  if (Front > 1 & W[Front] > 3*tiny)  					W[Front -1] <- tiny
  if (End < length(V))   												W[End +1] <- W[End]/2
  if ((End+1) < length(V) & W[End+1] > 3*tiny)  W[End+2] <- tiny
  for (x in 1:3) {  # in case neighbouring lowMid's
    for (i in seq_along(V)) {
      if (lowMid[i]) {
        W[i] <- (W[i-1] + W[i+1])/2
      }
    }
  }
  return( W )
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title LLR-age from Ageprior Matrix
#'
#' @description Get log10-likelihood ratios for a specific age difference from
#'   matrix \code{AgePriorExtra}.
#'
#' @param AgePriorExtra matrix in \code{\link{sequoia}} output
#' @param agedif vector with age differences, in whole numbers. Must occur in
#'   rownames of \code{AgePriorExtra}.
#' @param patmat numeric vector; choose maternal (1), paternal (2) relatives, or
#'   for each relationship the most-likely alternative (3).
#'
#' @return A matrix with \code{nrow} equal to the length of \code{agedif}, and 7
#'   columns: PO-FS-HS-GP-FA-HA-U.
#'
#' @details  This is a simple helper function to extract values from
#'  \code{AgePriorExtra}, e.g. to use together with \code{\link{CalcPairLL}}.
#'
#' @examples
#' # For a pair with unknown age difference, explore the difference age-based
#' # LLRs for all relationships, for a range of plausible age differences.
#' PairsG <- data.frame(ID1 = 'A', ID2 = 'B', AgeDif = rep(c(-2,2,3),2),
#'                    PatMat = rep(1:2, each=3))
#' cbind(PairsG,
#'       GetLLRAge(SeqOUT_griffin$AgePriorExtra,
#'                 agedif = PairsG$AgeDif, patmat = PairsG$PatMat))
#'
#' @export

GetLLRAge <- function(AgePriorExtra, agedif, patmat) {
  a <- as.character(agedif)
  k <- patmat
  if (!all(k %in% c(1,2,3)))  stop("patmat must be 1 (mat), 2 (pat), or 3 (unk)")
  if (!all(a %in% rownames(AgePriorExtra)))  stop("some agedif not in AgePriorExtra")
  if (length(a) != length(k))  stop("agedif and patmat must have same length")

  R <- list("M" = list("M", "FS", "MS", c("MGM", "MGF"), "MFA", c("MMA", "MPA")),
            "P" = list("P", "FS", "PS", c("PGM", "PGF"), "PFA", c("PMA", "PPA")))

  ALR <- matrix(0, length(a), 7,
                dimnames = list(seq_along(a),
                                c("PO", "FS", "HS", "GP", "FA", "HA", "U")))
  for (i in seq_along(a)) {
    ALRtmp <- matrix(NA, 2, 6,
                     dimnames=list(c("mat", "pat"),
                                   c("PO", "FS", "HS", "GP", "FA", "HA")))
    for (x in 1:2) {
      if (k[i]<3 & x!=k[i])  next
      ALRtmp[x,] <- log10(sapply(1:6, function(r) mean(AgePriorExtra[a[i], R[[x]][[r]]])))
    }
    if (k[i]==3) {
      ALR[i, 1:6] <- apply(ALRtmp, 2, max)
    } else {
      ALR[i, 1:6] <- ALRtmp[k[i],]
    }
  }

  return( round(ALR, 2) )
}
