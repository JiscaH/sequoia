#' @title Pedigree Reconstruction
#'
#' @description Perform pedigree reconstruction based on SNP data, including
#'   parentage assignment and sibship clustering.
#'
#' @details For each pair of candidate relatives, the likelihoods are calculated
#'   of them being parent-offspring (PO), full siblings (FS), half siblings
#'   (HS), grandparent-grandoffspring (GG), full avuncular (niece/nephew -
#'   aunt/uncle; FA), half avuncular/great-grandparental/cousins (HA), or
#'   unrelated (U). Assignments are made if the likelihood ratio (LLR) between
#'   the focal relationship and the most likely alternative exceed the threshold
#'   Tassign.
#'
#'   Dummy parents of sibships are denoted by F0001, F0002, ... (mothers)
#'   and M0001, M0002, ... (fathers), are appended to the bottom of the
#'   pedigree, and may have been assigned real or dummy parents themselves (i.e.
#'   sibship-grandparents). A dummy parent is not assigned to singletons.
#'
#'   Full explanation of the various options and interpretation of the output is
#'   provided in the vignettes and on the package website,
#'   https://jiscah.github.io/index.html .
#'
#' @param GenoM  numeric matrix with genotype data: One row per individual,
#'   one column per SNP, coded as 0, 1, 2, missing values as a negative number
#'   or NA. You can reformat data with \code{\link{GenoConvert}}, or use other
#'   packages to get it into a genlight object and then use \code{as.matrix}.
#' @param LifeHistData data.frame with up to 6 columns:
#'  \describe{
#'  \item{ID}{max. 30 characters long}
#'  \item{Sex}{1 = female, 2 = male, 3 = unknown, 4 = hermaphrodite,
#'            other numbers or NA = unknown}
#' \item{BirthYear }{birth or hatching year, integer, with missing values as NA
#'   or any negative number.}
#' \item{BY.min}{minimum birth year, only used if BirthYear is missing}
#' \item{BY.max}{maximum birth year, only used if BirthYear is missing}
#' \item{Year.last}{Last year in which individual could have had offspring. Can
#'   e.g. in mammals be the year before death for females, and year after death
#'   for males. } }
#' "Birth year" may be in any arbitrary discrete time unit relevant to the
#' species (day, month, decade), as long as parents are never born in the same
#' time unit as their offspring, and only integers are used. Individuals do not
#' need to be in the same order as in `GenoM', nor do all genotyped individuals
#' need to be included.
#' @param SeqList list with output from a previous run, to be re-used in the
#'   current run. Used are elements `PedigreePar', `LifeHist', `AgePriors',
#'   `Specs', and `ErrM', and these override the corresponding input parameters.
#'   Not all of these elements need to be present, and all other elements are
#'   ignored. If \code{SeqList$Specs} is provided, all  input parameters with
#'   the same name as its items are ignored, except \code{Module}.
#' @param Module one of
#'   \describe{
#'     \item{pre}{Only input check, return \code{SeqList$Specs}}
#'     \item{dup}{Also check for duplicate genotypes}
#'     \item{par}{Also perform parentage assignment (genotyped parents to
#'       genotyped offspring)}
#'     \item{ped}{(Also) perform full pedigree reconstruction, including
#'       sibship clustering and grandparent assignment. By far the most time
#'       consuming, and may take several hours for large datasets.}
#'   }
#' @param Err assumed per-locus genotyping error rate, as a single number, or a
#'   length 3 vector with P(hom|hom), P(het|hom), P(hom|het), or a 3x3 matrix.
#'   See details below. The error rate is presumed constant across SNPs, and
#'   missingness is presumed random with respect to actual genotype. Using
#'   \code{Err} >5\% is not recommended, and \code{Err} >10\% strongly
#'   discouraged. See \code{\link{Err_RADseq}} to convert per-allele rates at
#'   homozygous and heterozygous sites to the required length-3 vector, and
#'   \code{\link{ErrToM}} for further genotyping error details.
#' @param ErrFlavour function that takes \code{Err} (single number) as input,
#'   and returns a length 3 vector or 3x3 matrix, or choose from inbuilt options
#'   'version2.9', 'version2.0', 'version1.3', or 'version1.1', referring to the
#'   sequoia version in which they were the default. Ignored if \code{Err} is a
#'   vector or matrix.
#' @param Tfilter threshold log10-likelihood ratio (LLR) between a proposed
#'   relationship versus unrelated, to select candidate relatives. Typically a
#'   negative value, related to the fact that unconditional likelihoods are
#'   calculated during the filtering steps. More negative values may decrease
#'   non-assignment, but will increase computational time.
#' @param Tassign minimum LLR required for acceptance of proposed relationship,
#'   relative to next most likely relationship. Higher values result in more
#'   conservative assignments. Must be zero or positive.
#' @param MaxSibshipSize  maximum number of offspring for a single individual (a
#'   generous safety margin is advised).
#' @param DummyPrefix character vector of length 2 with prefixes for dummy dams
#'   (mothers) and sires (fathers); maximum 20 characters each. Length 3 vector
#'   in case of hermaphrodites (or default prefix 'H').
#' @param Complex  Breeding system complexity. Either "full" (default), "simp"
#'   (simplified, no explicit consideration of inbred relationships), "mono"
#'   (monogamous).
#' @param Herm  Hermaphrodites, either "no", "A" (distinguish between dam and
#'   sire role, default if at least 1 individual with sex=4), or "B" (no
#'   distinction between dam and sire role). Both of the latter deal with
#'   selfing.
#' @param UseAge  either "yes" (default), "no" (only use age differences for
#'   filtering), or "extra" (additional rounds with extra reliance on ageprior,
#'   may boost assignments but increased risk of erroneous assignments). Used
#'   during full reconstruction only.
#' @param args.AP list with arguments to be passed on to
#'   \code{\link{MakeAgePrior}}, e.g. `Discrete` (non-overlapping generations),
#'   `MinAgeParent`, `MaxAgeParent`.
#' @param mtSame  \strong{NEW} matrix indicating whether individuals (might)
#'   have the same mitochondrial haplotype (1), and may thus be matrilineal
#'   relatives, or not (0). Row names and column names should match IDs in
#'   `GenoM`. Not all individuals need to be included and order is not
#'   important. Please report any issues. For details see the mtDNA vignette.
#' @param CalcLLR  TRUE/FALSE; calculate log-likelihood ratios for all assigned
#'   parents (genotyped + dummy; parent vs. otherwise related). Time-consuming
#'   in large datasets. Can be done separately with \code{\link{CalcOHLLR}}.
#' @param quiet suppress messages: TRUE/FALSE/"verbose".
#' @param Plot display plots from \code{\link{SnpStats}, \link{MakeAgePrior}},
#'   and \code{\link{SummarySeq}}. Defaults (NULL) to TRUE when quiet=FALSE or
#'   "verbose", and FALSE when quiet=TRUE. If you get error 'figure margins too
#'   large', enlarge the plotting area (drag with mouse). Error 'invalid
#'   graphics state' can be dealt with by clearing the plotting area with
#'   dev.off().
#' @param StrictGenoCheck Automatically exclude any individuals genotyped for
#'   <5% of SNPs, and any SNPs genotyped for <5% of individuals (TRUE); this was
#'   the unavoidable default up to version 2.4.1. Otherwise only excluded are
#'   (very nearly) monomorphic SNPs, SNPs scored for fewer than 2 individuals,
#'   and individuals scored for fewer than 2 SNPs.
#'
#'
#' @return A list with some or all of the following components, depending on
#'   \code{Module}. All input except \code{GenoM} is included in the output.
#' \item{AgePriors}{Matrix with age-difference based probability ratios for
#'   each relationship, used for full pedigree reconstruction; see
#'   \code{\link{MakeAgePrior}} for details. When running only parentage
#'   assignment (\code{Module="par"}) the returned AgePriors has been updated to
#'   incorporate the information of the assigned parents, and is ready for use
#'   during full pedigree reconstruction.}
#' \item{args.AP}{(input) arguments used to specify age prior matrix. If a
#'   custom ageprior was provided via \code{SeqList$AgePrior}, this matrix is
#'   returned instead}
#' \item{DummyIDs}{Dataframe with pedigree for dummy individuals, as well as
#' their sex, estimated birth year (point estimate, upper and lower bound of
#' 95\% confidence interval; see also \code{\link{CalcBYprobs}}), number of
#' offspring, and offspring IDs. From version 2.1 onwards, this includes dummy
#' offspring.}
#' \item{DupGenotype}{Dataframe, duplicated genotypes (with different IDs,
#'  duplicate IDs are not allowed). The specified number of maximum mismatches
#'   is used here too. Note that this dataframe may include pairs of closely
#'   related individuals, and monozygotic twins.}
#' \item{DupLifeHistID}{Dataframe, row numbers of duplicated IDs in life
#'   history dataframe. For convenience only, but may signal a problem. The
#'   first entry is used.}
#' \item{ErrM}{(input) Error matrix; probability of observed genotype (columns)
#'   conditional on actual genotype (rows)}
#' \item{ExcludedInd}{Individuals in GenoM which were excluded because of a
#'   too low genotyping success rate (<50\%).}
#' \item{ExcludedSNPs}{Column numbers of SNPs in GenoM which were excluded
#'   because of a too low genotyping success rate (<10\%).}
#' \item{LifeHist}{(input) Dataframe with sex and birth year data. All missing
#'   birth years are coded as '-999', all missing sex as '3'.}
#' \item{LifeHistPar}{LifeHist with additional columns 'Sexx' (inferred Sex when
#' assigned as part of parent-pair), 'BY.est' (mode of birth year probability
#' distribution), 'BY.lo' (lower limit of 95\% highest density region), 'BY.hi'
#' (higher limit), inferred after parentage assignment. 'BY.est' is NA when the
#' probability distribution is flat between 'BY.lo' and 'BY.hi'.}
#' \item{LifeHistSib}{as LifeHistPar, but estimated after full pedigree
#' reconstruction}
#' \item{NoLH}{Vector, IDs in genotype data for which no life history data is
#'  provided.}
#' \item{Pedigree}{Dataframe with assigned genotyped and dummy parents from
#'   Sibship step; entries for dummy individuals are added at the bottom.}
#' \item{PedigreePar}{Dataframe with assigned parents from Parentage step.}
#' \item{Specs}{Named vector with parameter values.}
#' \item{TotLikParents}{Numeric vector, Total likelihood of the genotype data
#'   at initiation and after each iteration during Parentage.}
#' \item{TotLikSib}{Numeric vector, Total likelihood of the genotype data
#'   at initiation and after each iteration during Sibship clustering.}
#' \item{AgePriorExtra}{As AgePriors, but including columns for grandparents
#'  and avuncular pairs. NOT updated after parentage assignment, but returned
#'  as used during the run.}
#' \item{DummyClones}{Hermaphrodites only: female-male dummy ID pairs that refer
#'   to the same non-genotyped individual}
#'
#' List elements PedigreePar and Pedigree both have the following columns:
#'  \item{id}{Individual ID}
#'  \item{dam}{Assigned mother, or NA}
#'  \item{sire}{Assigned father, or NA}
#'  \item{LLRdam}{Log10-Likelihood Ratio (LLR) of this female being the mother,
#'  versus the next most likely relationship between the focal individual and
#'  this female. See Details below for relationships considered, and see
#'  \code{\link{CalcPairLL}} for underlying likelihood values and further
#'  details)}
#'  \item{LLRsire}{idem, for male parent}
#'  \item{LLRpair}{LLR for the parental pair, versus the next most likely
#'   configuration between the three individuals (with one or neither parent
#'   assigned)}
#'  \item{OHdam}{Number of loci at which the offspring and mother are
#'    opposite homozygotes}
#'  \item{OHsire}{idem, for father}
#'  \item{MEpair}{Number of Mendelian errors between the offspring and the
#'    parent pair, includes OH as well as e.g. parents being opposing
#'    homozygotes, but the offspring not being a heterozygote. The offspring
#'    being OH with both parents is counted as 2 errors.}
#'
#'
#' @section Genotyping error rate:
#'   The genotyping error rate \code{Err} can be specified three different ways:
#'   \itemize{
#'    \item A single number, which is combined with \code{ErrFlavour} by
#'      \code{\link{ErrToM}} to create a length 3 vector (next item). By
#'      default (\code{ErrFlavour} = 'version2.9'), P(hom|hom)=$(E/2)^2$,
#'      P(het|hom)=$E-(E/2)^2$, P(hom|het)=$E/2$.
#'   \item a length 3 vector (NEW from version 2.6), with the probabilities to
#'     observe a actual homozygote as the other homozygote (hom|hom), to observe
#'     a homozygote as heterozygote (het|hom), and to observe an actual
#'     heterozygote as homozygote (hom|het). This assumes that the two alleles
#'     are equivalent with respect to genotyping errors, i.e. $P(AA|aa) =
#'     P(aa|AA)$, $P(aa|Aa)=P(AA|Aa)$, and $P(aA|aa)=P(aA|AA)$.
#'   \item a 3x3 matrix, with the probabilities of observed genotype (columns)
#'     conditional on actual genotype (rows). Only needed when the assumption
#'     in the previous item does not hold. See \code{\link{ErrToM}} for details.
#'   }
#'
#' @section (Too) Few Assignments?:
#' Possibly \code{Err} is much lower than the actual genotyping error rate.
#'
#' Alternatively, a true parent will not be assigned when it is:
#' \itemize{
#'   \item unclear who is the parent and who the offspring, due to unknown birth
#'   year for one or both individuals
#'   \item unclear whether the parent is the father or mother
#'   \item unclear if it is a parent or e.g. full sibling or grandparent, due to
#'   insufficient genetic data
#'   }
#'  And true half-siblings will not be clustered when it is:
#'  \itemize{
#'    \item unclear if they are maternal or paternal half-siblings
#'    \item unclear if they are half-siblings, full avuncular, or grand-parental
#'    \item unclear what type of relatives they are due to insufficient genetic
#'    data
#'  }
#'  All pairs of non-assigned but likely/definitely relatives can be found with
#'  \code{\link{GetMaybeRel}} (run without conditioning on a pedigree). For a
#'  method to do pairwise 'assignments', see
#'  https://jiscah.github.io/articles/pairLL_classification.html ; for further
#'  information, see the vignette.
#'
#'  If you already had a partial pedigree, running \code{\link{CalcOHLLR}} or
#'  \code{\link{CalcParentProbs}} on it with the new SNP data and various
#'  parameter value combinations may be informative.
#'
#' @section Disclaimer:
#' While every effort has been made to ensure that sequoia provides what it
#' claims to do, there is absolutely no guarantee that the results provided are
#' correct. Use of sequoia is entirely at your own risk.
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @references Huisman, J. (2017) Pedigree reconstruction from SNP data:
#'   Parentage assignment, sibship clustering, and beyond. Molecular Ecology
#'   Resources 17:1009--1024.
#'
#' @section Website:
#' https://jiscah.github.io/
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{GenoConvert}} to read in various data formats,
#'   \item \code{\link{CheckGeno}}, \code{\link{SnpStats}} to calculate
#'     missingness and allele frequencies,
#'   \item \code{\link{SimGeno}}  to simulate SNP data from a pedigree,
#'   \item \code{\link{MakeAgePrior}} to estimate effect of age on relationships,
#'   \item \code{\link{GetMaybeRel}} to find pairs of potential relatives,
#'   \item \code{\link{SummarySeq}} and \code{\link{PlotAgePrior}} to visualise
#'   results,
#'   \item \code{\link{GetRelM}} to turn a pedigree into pairwise relationships,
#'   \item \code{\link{CalcOHLLR}} to calculate Mendelian errors and LLR for any
#'    pedigree,
#'   \item \code{\link{CalcParentProbs}} to calculate assignment probabilities
#'   instead of LLRs,
#'   \item \code{\link{CalcPairLL}} for likelihoods of various relationships
#'    between specific pairs,
#'   \item \code{\link{CalcBYprobs}} to estimate birth years,
#'   \item \code{\link{PedCompare}} and \code{\link{ComparePairs}} to compare to
#'   two pedigrees,
#'   \item \code{\link{EstConf}} to estimate assignment errors,
#'   \item \code{\link{writeSeq}} to save results,
#'   \item \code{vignette("sequoia")} for detailed manual & FAQ.
#' }
#'
#'
#'
#' @examples
#' # ===  EXAMPLE 1: simulated data  ===
#' head(SimGeno_example[,1:10])
#' head(LH_HSg5)
#' # parentage assignment:
#' SeqOUT <- sequoia(GenoM = SimGeno_example, Err = 0.005,
#'                   LifeHistData = LH_HSg5, Module="par", Plot=TRUE)
#' names(SeqOUT)
#' SeqOUT$PedigreePar[34:42, ]
#'
#' # compare to true (or old) pedigree:
#' PC <- PedCompare(Ped_HSg5, SeqOUT$PedigreePar)
#' PC$Counts["GG",,]
#'
#' \donttest{
#' # parentage assignment + full pedigree reconstruction:
#' # (note: this can be rather time consuming)
#' SeqOUT2 <- sequoia(GenoM = SimGeno_example, Err = 0.005,
#'                   LifeHistData = LH_HSg5, Module="ped", quiet="verbose")
#' SeqOUT2$Pedigree[34:42, ]
#'
#' PC2 <- PedCompare(Ped_HSg5, SeqOUT2$Pedigree)
#' PC2$Counts["GT",,]
#' PC2$Counts[,,"dam"]
#'
#' # different kind of pedigree comparison:
#' ComparePairs(Ped1=Ped_HSg5, Ped2=SeqOUT$PedigreePar, patmat=TRUE)
#'
#' # results overview:
#' SummarySeq(SeqOUT2)
#'
#' # important to run with approx. correct genotyping error rate:
#' SeqOUT2.b <- sequoia(GenoM = SimGeno_example, #  Err = 1e-4 by default
#'                   LifeHistData = LH_HSg5, Module="ped", Plot=FALSE)
#' PC2.b <- PedCompare(Ped_HSg5, SeqOUT2.b$Pedigree)
#' PC2.b$Counts["GT",,]
#' }
#'
#' \dontrun{
#' # ===  EXAMPLE 2: real data  ===
#' # ideally, select 400-700 SNPs: high MAF & low LD
#' # save in 0/1/2/NA format (PLINK's --recodeA)
#' GenoM <- GenoConvert(InFile = "inputfile_for_sequoia.raw",
#'                      InFormat = "raw")  # can also do Colony format
#' SNPSTATS <- SnpStats(GenoM)
#' # perhaps after some data-cleaning:
#' write.table(GenoM, file="MyGenoData.txt", row.names=T, col.names=F)
#'
#' # later:
#' GenoM <- as.matrix(read.table("MyGenoData.txt", row.names=1, header=FALSE))
#' # or for very large datasets:
#' GenoM <- data.table::fread('MyGenoData.txt') %>% as.matrix(rownames=1)
#'
#' LHdata <- read.table("LifeHistoryData.txt", header=T) # ID-Sex-birthyear
#' SeqOUT <- sequoia(GenoM, LHdata, Err=0.005)
#' SummarySeq(SeqOUT)
#'
#' SeqOUT$notes <- "Trial run on cleaned data"  # add notes for future reference
#' saveRDS(SeqOUT, file="sequoia_output_42.RDS")  # save to R-specific file
#' writeSeq(SeqOUT, folder="sequoia_output")  # save to several plain text files
#'
#' # runtime:
#' SeqOUT$Specs$TimeEnd - SeqOUT$Specs$TimeStart
#' }
#'
#' @export

sequoia <- function(GenoM = NULL,
                    LifeHistData = NULL,
                    SeqList = NULL,
                    Module = "ped",
                    Err = 0.0001,
                    Tfilter = -2.0,
                    Tassign = 0.5,
                    MaxSibshipSize = 100,
                    DummyPrefix = c("F", "M"),
                    Complex = "full",
                    Herm = "no",
                    UseAge = "yes",
                    args.AP = list(Flatten=NULL, Smooth=TRUE),
                    mtSame = NULL,
                    CalcLLR = FALSE,
                    quiet = FALSE,
                    Plot = NULL,
                    StrictGenoCheck = TRUE,
                    ErrFlavour = "version2.9")
{
  TimeStart <- Sys.time()

  # set quiet & Plot ----
  if (!quiet %in% c(TRUE, FALSE, "verbose"))
    stop("'quiet' must be TRUE or FALSE or 'verbose'")
  quietR <- ifelse(quiet == "verbose", FALSE, quiet)  # 'verbose' determines chattiness of Fortran only
  if (is.null(Plot))   # default
    Plot <- ifelse(quietR, FALSE, TRUE)

  # check input ----
  if (!is.null(LifeHistData) & !inherits(LifeHistData, 'data.frame'))
    stop("`LifeHistData` must be a data.frame or NULL")
  if (!is.null(SeqList) & !inherits(SeqList, 'list'))
    stop("`SeqList` must be a list or `NULL`")

  if (!is.null(SeqList)) {
    SeqList_names <- c("Specs", "ErrM", "args.AP", "AgePriors", "LifeHist","PedigreePar")
    if (!any(names(SeqList) %in% SeqList_names) | any(is.na(names(SeqList))) ) {
      cli::cli_alert_danger("You seem to have misspelled one or more names of elements of `SeqList`:")
      cli::cli_li(setdiff(names(SeqList), SeqList_names))
      stop("Please correct `SeqList` names")
    }
  }


  # Module ----
  # MaxSibIter is now fully deprecated.
  Module <- factor(Module, levels = c("pre", "dup", "par", "ped"))
  if (is.na(Module))  stop("'Module'  must be 'pre', 'dup', 'par', or 'ped'")
  dostep <- c('pre' = TRUE,
              'dup' = Module!='pre',
              'par' = Module %in% c('par','ped'),
              'ped' = Module=='ped')
  if ('DupGenotype' %in% names(SeqList))  dostep['dup'] <- FALSE
  if ('PedigreePar' %in% names(SeqList) & Module=='ped')  dostep['par'] <- FALSE
  # module='par' + PedigreePar: parentage assignment with pedigree-prior
  # for output:
  DupList <- NULL
  ParList <- NULL
  SibList <- NULL

  # Check genotype matrix ----
  GenoM[is.na(GenoM)] <- -9
  if (!quietR) cli::cli_alert_info("Checking input data ...")
  Excl <- CheckGeno(GenoM, quiet=quietR, Plot=Plot, Return = "excl",
                    Strict = StrictGenoCheck, DumPrefix=DummyPrefix)
  if ("ExcludedSnps" %in% names(Excl))  GenoM <- GenoM[, -Excl[["ExcludedSnps"]]]
  if ("ExcludedSnps-mono" %in% names(Excl))  GenoM <- GenoM[, -Excl[["ExcludedSnps-mono"]]]
  if ("ExcludedIndiv" %in% names(Excl))  GenoM <- GenoM[!rownames(GenoM) %in% Excl[["ExcludedIndiv"]], ]


  # Check life history data ----
  if ("LifeHist" %in% names(SeqList)) {
    if (!quietR)  cli::cli_alert_info("using `LifeHistData` in `SeqList`")
    LifeHistData <- SeqList$LifeHist
  } else if (is.null(LifeHistData) & Module != "dup") {
    cli::cli_alert_warning("no `LifeHistData` provided, expect lower assignment rate")
  }
  ChkLH.L <- CheckLH(LifeHistData, gID = rownames(GenoM),  sorted=FALSE,
                     returnDups = TRUE)
  # keep non-genotyped IDs (for future reference); orderLH() called by SeqParSib()
  OUT_LH <- ChkLH.L[c("DupLifeHistID", "NoLH")]
  LifeHistData <- ChkLH.L$LifeHistData   # duplicates removed, if any

  if (!quietR) {
    cli::cli_h3("Among genotyped individuals: ___")
    gID <- rownames(GenoM)
    tbl_sex <- table(factor(LifeHistData$Sex[LifeHistData$ID %in% gID], levels=1:4))
    cli::cli_alert_info(c("There are {tbl_sex['1']} females, {tbl_sex['2']} males, ",
            "{tbl_sex['3']} of unknown sex, and {tbl_sex['4']} hermaphrodites."))
    range_Year <- matrix(NA,4,2, dimnames=list(c("BirthYear", "BY.min", "BY.max", 'Year.last'),
                                            c('min', 'max')))
    for (x in rownames(range_Year)) {
      range_Year[x,] <- suppressWarnings(
        range(LifeHistData[,x][LifeHistData$ID %in% gID & LifeHistData[,x] >= 0]))
    }
    if (any(is.finite(range_Year[1,])))
      cli::cli_alert_info("Exact birth years are from {range_Year[1,1]} to {range_Year[1,2]}")
    if (any(is.finite(range_Year[2:3,])))
      cli::cli_alert_info(c("Birth year min/max are from {min(range_Year[2:3,1], na.rm=TRUE)} to ",
            "{max(range_Year[2:3,2], na.rm=TRUE)}"))
    if (all(!is.finite(range_Year)))
      cli::cli_alert_info("All birth years are unknown")
    if (any(is.finite(range_Year[4,])))
      cli::cli_alert_info("`Year.Last` are from {range_Year[4,1]} to {range_Year[4,2]}")
    cli::cli_text("___\n")
  }

  utils::flush.console()    # print all warnings thus far

  # Check/make ageprior ----
  if ("AgePriors" %in% names(SeqList)) {
    if(!quietR)  cli::cli_alert_info("using `AgePriors` in `SeqList`")
    AgePriors <- CheckAP( SeqList[["AgePriors"]] )
  } else {
    if(!quietR)  cli::cli_alert_info("Calling `MakeAgePrior()` ...")
    AgePriors <- do.call(MakeAgePrior, c(list(Pedigree = SeqList[["PedigreePar"]],
                    LifeHistData = LifeHistData[LifeHistData$ID %in% rownames(GenoM),],
                                              Plot = Plot,
                                              quiet = ifelse(Module=="dup", TRUE, quietR)),
                                         args.AP))
  }

  utils::flush.console()

  # check & reformat mitochondrial same/not matrix
  mtDif <- mtSame2Dif(mtSame, gID = rownames(GenoM))

  # make list with parameter values (loose input or SeqList$Specs) ----
  if ("Specs" %in% names(SeqList)) {
    PARAM <- SpecsToParam(SeqList$Specs, SeqList$ErrM, ErrFlavour,  # re-wrapping + ErrToM check
                          dimGeno = dim(GenoM), Module, quiet) # overrule old values
  } else {
    if (is.logical(UseAge))   UseAge <- ifelse(UseAge, "yes", "no")
    if ((Herm != "no" | any(LifeHistData$Sex==4, na.rm=TRUE)) & length(DummyPrefix)==2)
      DummyPrefix <- c(DummyPrefix, "H")

    PARAM <- namedlist(dimGeno = dim(GenoM),
                       Err,
                       ErrFlavour,
                       Tfilter,
                       Tassign,
                       nAgeClasses = nrow(AgePriors),
                       MaxSibshipSize,
                       Module = as.character(Module),
                       DummyPrefix,
                       Complex,
                       Herm,
                       UseAge,
                       CalcLLR,
                       quiet)
    PARAM$ErrM <- ErrToM(Err, flavour = ErrFlavour, Return = "matrix")
  }

  # MaxMismatch ----
  # vector with max. mismatches for duplicates, PO pairs, PPO trios
  if (!"MaxMismatchV" %in% names(PARAM)) {  # DUP/OH/ME from version 2.0 onwards
    sts <- SnpStats(GenoM, Plot=FALSE)
    PARAM$MaxMismatchV <- setNames(CalcMaxMismatch(Err=PARAM$ErrM,
                                                   MAF=sts[,"AF"],
                                                   ErrFlavour=PARAM$ErrFlavour,
                                                   qntl=0.9999^(1/nrow(GenoM))),
                                   c("DUP", "OH", "ME"))
  }


  # hermaphrodites ----
  if (any(LifeHistData$Sex==4, na.rm=TRUE) && PARAM$Herm == "no") {  #!grepl("herm", PARAM$Complex)) {
    if (!quietR) cli::cli_alert_warning("Detected hermaphrodites (sex=4), changing `Herm` to 'A'")
#    PARAM$Complex <- "herm"
    PARAM$Herm <- "A"
  }
  if (PARAM$Herm == "B" && any(LifeHistData$Sex %in% c(1,2)))
     cli::cli_alert_warning("Results may be inconsistent when combining `Herm='B'` with known-sex individuals")


  # check parameter values ----
  CheckParams(c(PARAM, list(Plot=Plot)))
  utils::flush.console()

  # turn into fortran-friendly parameters ----
  FortPARAM <- MkFortParams(PARAM, fun="main")  # SpecsInt, SpecsDbl, character choice to numeric


  # @@ 2 @@ Duplicate check ----
  if (dostep[['dup']]) {
    if(!quietR)  message(cli::col_green("\n~~~ Duplicate check ~~~"))
    DupList <- DuplicateCheck(GenoM, FortPARAM, quiet=quietR)
    utils::flush.console()
    if ("DupGenoID" %in% names(DupList)) {  # fatal error  (message by DuplicateCheck())
      return(DupList)
    } else if (length(DupList)==0 & !quietR) {
      cli::cli_alert_success("No potential duplicates found")
    } else if ('DupGenotype' %in% names(DupList) && nrow(DupList[['DupGenotype']])>=nrow(GenoM)) {
      Sys.sleep(1)   # ensures warnings are printed first
      cli::cli_alert_danger(cli::col_red(paste('The high number of putative duplicate genotypes suggests',
                                  'that the genotype data is not informative enough for',
                                  'parentage assignment. Interpret results with caution!')), wrap=TRUE)
      # dostep[c('par', 'ped')] <- FALSE
    }
  }
  utils::flush.console()


  # @@ 3 @@ Parentage assignment ----
  if (any(dostep[c('par', 'ped')]) & "PedigreePar" %in% names(SeqList)) {
    if (!quietR) {
      message(cli::col_green("\n~~~ Parentage assignment ~~~"))
      cli::cli_alert_info("using `PedigreePar` in `SeqList`")
    }
    PedParents <- PedPolish(SeqList[["PedigreePar"]], gID = rownames(GenoM),
                            ZeroToNA = TRUE, DropNonSNPd = TRUE)
  } else {
    PedParents <- NULL
  }

  if (dostep[['par']]) {
    if(!quietR & !is.null(PedParents)) {
      message(cli::col_green("\n~~~ Parentage assignment with pedigree-prior ~~~"))  # only sensible with Herm=A or B (?)
    } else if(!quietR) {
      message(cli::col_green("\n~~~ Parentage assignment ~~~"))
    }

    ParList <- SeqParSib(ParSib = "par", FortPARAM, GenoM,
                         LhIN=LifeHistData, AgePriors=AgePriors,
                         Parents=PedParents, mtDif=mtDif,
                         DumPfx = PARAM$DummyPrefix, quiet=quietR)
    if (Plot) {
      SummarySeq(ParList, Panels="G.parents")
    }

  } else if (dostep[['ped']]) {  # PedigreePar provided as input
    ParList <- list(PedigreePar = PedParents)   # re-use assigned parents given as SeqList input (after polishing)
  }
  utils::flush.console()

  # check that PedigreePar is a valid pedigree (no indiv is its own ancestor)
  W <- tryCatch.W.E(getGenerations(ParList$PedigreePar, StopIfInvalid=FALSE))$warning
  if (!is.null(W) & dostep[['ped']]) {
    dostep['ped'] <- FALSE
  }


  # Update ageprior ----
  if (any(dostep[c('par', 'ped')]) & !"AgePriors" %in% names(SeqList) & is.null(W) &
      !"PedigreePar" %in% names(SeqList) && !is.null(LifeHistData)) {
    if (Plot)  Sys.sleep(1)  # else no time to notice previous plot
    # in case MakeAgePrior throws error, do return parentage results:
    AgePriors <- tryCatch( do.call(MakeAgePrior,
                                   c(list(Pedigree = ParList$PedigreePar[, 1:3],
                                          LifeHistData = LifeHistData[LifeHistData$ID %in% rownames(GenoM),],
                                          Plot = Plot & dostep['ped'],
                                          quiet = !(!quietR & dostep['ped'])),
                                     args.AP)),
                           error = function(e) {
                             message("AgePrior error! \n", e)
                             return(NA)
                           })
    if (all(is.na(AgePriors)))  return(ParList)

  } else if ("AgePriors" %in% names(SeqList) & !"PedigreePar" %in% names(SeqList)) {
    if(!quietR)  cli::cli_alert_info("using `AgePriors` in `SeqList` again")
  }

  if (nrow(AgePriors) != PARAM$nAgeClasses) {
    PARAM$nAgeClasses <- nrow(AgePriors)
    FortPARAM$SpecsInt[["nAgeCl"]] <- nrow(AgePriors)
  }


  # @@ 4 @@ Full pedigree reconstruction ----
  if (dostep['ped']) {
    if (!all(apply(AgePriors, 2, function(x) any(x > 0))))
      stop("AgePriors error: some relationships are impossible for all age differences")
    if(!quietR)  message(cli::col_green("\n~~~ Full pedigree reconstruction ~~~"))
    SibList <- SeqParSib(ParSib = "sib", FortPARAM, GenoM,
                         LhIN = LifeHistData, AgePriors = AgePriors,
                         Parents = ParList$PedigreePar, mtDif=mtDif,
                         DumPfx = PARAM$DummyPrefix, quiet = quietR)
    ParList <- ParList[names(ParList) != "AgePriorExtra"]  # else included 2x w same name
  }

  #=====================
  # Output ----

  if (quiet=='verbose' & any(dostep[c('par', 'ped')])) {
    cli::cli_alert_info(paste("You can use `SummarySeq()` for pedigree details, ",
                        "`CalcParentProbs()` for assignment probabilities, ",
                        "and `EstConf()` for confidence estimates"), wrap=TRUE)
    cli::cli_alert_info(paste('Run `GetMaybeRel()` conditional on this pedigree to check for any ',
                              'non-assigned ', c(par = 'parents', ped = 'relatives')[as.character(Module)]),
                              wrap=TRUE)
  }

  OUT <- list()
  OUT[["Specs"]] <- ParamToSpecs(PARAM, TimeStart, ErrFlavour)
  OUT[["ErrM"]] <- PARAM$ErrM
  if (is.function(ErrFlavour))  OUT[["ErrFlavour"]] <- ErrFlavour
  if ("AgePriors" %in% names(SeqList)) {
    OUT[["args.AP"]] <- SeqList[['AgePriors']]
  } else {
    OUT[["args.AP"]] <- args.AP
  }
  if (length(Excl)>0)  OUT <- c(OUT, Excl)
  if (length(OUT_LH)>0)  OUT <- c(OUT, OUT_LH)
  OUT <- c(OUT,
           list(AgePriors = AgePriors,
                LifeHist = LifeHistData),
           DupList, ParList, SibList)

  if (Plot & dostep[['ped']]) {
    PlotAgePrior(OUT$AgePriorExtra)
    SummarySeq(OUT, Panels="G.parents")
  }

  return(OUT)
}
