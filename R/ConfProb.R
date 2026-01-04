#============================================================================
#============================================================================
#' @title Confidence Probabilities
#'
#' @description Estimate confidence probabilities ('backward') and assignment
#'   error rates ('forward') per category (genotyped/dummy) by repeatedly
#'   simulating genotype data from a reference pedigree using
#'   \code{\link{SimGeno}}, reconstruction a pedigree from this using
#'   \code{\link{sequoia}}, and counting the number of mismatches using
#'   \code{\link{PedCompare}}.
#'
#' @details The confidence probability is taken as the number of correct
#'   (matching) assignments, divided by all assignments made in the
#'   \emph{observed} (inferred-from-simulated) pedigree. In contrast, the false
#'   negative & false positive assignment rates are proportions of the number of
#'   parents in the \emph{true} (reference) pedigree. Each rate is calculated
#'   separately for dams & sires, and separately for each category
#'   (\strong{G}enotyped/\strong{D}ummy(fiable)/\strong{X} (none)) of
#'   individual, parent and co-parent.
#'
#'  This function does not know which individuals in the actual \code{Pedigree}
#'  are genotyped, so the confidence probabilities need to be added to the
#'  \code{Pedigree} as shown in the example at the bottom.
#'
#'  A confidence of \eqn{1} means all assignments on simulated data were correct for
#'  that category-combination. It should be interpreted as (and perhaps modified
#'  to) \eqn{> 1 - 1/N}, where sample size \code{N} is given in the last column
#'  of the \code{ConfProb} and \code{PedErrors} dataframes in the output. The
#'  same applies for a false negative/positive rate of \eqn{0} (i.e. to be
#'  interpreted as \eqn{< 1/N}).
#'
#' @section Assumptions:
#'   Because the actual true pedigree is (typically) unknown, the provided
#'   reference pedigree is used as a stand-in and assumed to be the true
#'   pedigree, with unrelated founders. It is also assumed that the probability
#'   to be genotyped is equal for all parents; in each round, a new random
#'   set of parents (proportion set by \code{ParMis}) is mimicked to be
#'   non-genotyped. In addition, SNPs are assumed to segregate independently.
#'
#'   An experimental version offering more fine-grained control is available at
#'   https://github.com/JiscaH/sequoiaExtra .
#'
#' @section Object size:
#'   The size in Kb of the returned list can become pretty big, as each of the
#'   inferred pedigrees is included. When running \code{EstConf} many times for
#'   a range of parameter values, it may be prudent to save the required summary
#'   statistics for each run rather than the full output.
#'
#' @section Errors:
#'   If you have a large pedigree and try to run this function on multiple
#'   cores, you may run into "Cannot allocate vector of size ..." errors or even
#'   unexpected crashes: there is not enough computer memory for each separate
#'   run. Try reducing `nCores`.
#'
#' @param Pedigree reference pedigree from which to simulate, dataframe with
#'   columns id-dam-sire. Additional columns are ignored.
#' @param LifeHistData dataframe with id, sex (1=female, 2=male, 3=unknown),
#' birth year, and optionally BY.min - BY.max - YearLast.
#' @param args.sim  list of arguments to pass to \code{\link{SimGeno}}, such as
#'   \code{nSnp} (number of SNPs), \code{SnpError} (genotyping error rate) and
#'   \code{ParMis} (proportion of non-genotyped parents). Set to \code{NULL} to
#'   use all default values.
#' @param args.seq  list of arguments to pass to \code{\link{sequoia}}, such as
#'   \code{Module} ('par' or 'ped'), \code{Err} (assumed genotyping error rate),
#'   and \code{Complex}. May include (part of) \code{SeqList}, a list of sequoia
#'   output (i.e. as a list-within-a-list). Set to \code{NULL} to use all
#'   default values.
#' @param nSim number of rounds of simulate - reconstruct - compare to
#'   perform, i.e. number of simulated datasets.
#' @param nCores number of computer cores to use. If \code{>1}, package
#'   \pkg{parallel} is used. Set to NULL to use all but one of the available
#'   cores, as detected by \code{parallel::detectCores()} (using all cores tends
#'   to freeze up your computer). With large datasets, the amount of computer
#'   memory may be the limiting factor for the number of cores you can use.
#' @param quiet suppress messages. \code{TRUE} runs \code{SimGeno} and
#'   \code{sequoia} quietly, \code{'very'} also suppresses other messages and
#'   the simulation counter when \code{nCores=1} (there is no simulation counter
#'   when \code{nCores>1}).
#'
#' @return A list, with elements:
#'   \item{ConfProb}{See below}
#'   \item{PedErrors}{See below}
#'   \item{Pedigree.reference}{the pedigree from which data was simulated}
#'   \item{LifeHistData}{}
#'   \item{Pedigree.inferred}{a list with for each simulation the inferred
#'     pedigree based on the simulated data}
#'   \item{SimSNPd}{a list with for each simulation the IDs of the individuals
#'     simulated to have been genotyped}
#'   \item{PedComp.fwd}{array with \code{Counts} from the 'forward'
#'     \code{PedCompare}, from which \code{PedErrors} is calculated}
#'   \item{RunParams}{a list with the call to \code{EstConf} as a semi-nested
#'   list (args.sim, args.seq, nSim, nCores), as well as the default parameter
#'   values for \code{SimGeno} and \code{sequoia}.}
#'   \item{RunTime}{\code{sequoia} runtime per simulation in seconds, as
#'     measured by \code{\link{system.time}()['elapsed']}.}
#'
#' Dataframe \code{ConfProb} has 7 columns:
#' \item{id.cat, dam.cat, sire.cat}{Category of the focal individual, dam, and
#'   sire, in the pedigree inferred based on the simulated data. Coded as
#'   G=genotyped, D=dummy, X=none}
#' \item{dam.conf}{Probability that the dam is correct, given the categories of
#'   the assigned dam and sire (ignoring whether or not the sire is correct)}
#' \item{sire.conf}{as \code{dam.conf}, for the sire}
#' \item{pair.conf}{Probability that both dam and sire are correct, given their
#'   categories}
#' \item{N}{Number of individuals per category-combination, across all
#'   \code{nSim} simulations}
#'
#' Array \code{PedErrors} has three dimensions:
#' \item{class}{\itemize{
#'   \item \code{FalseNeg}(atives): could have been assigned but was not
#' (individual + parent both genotyped or dummifiable; P1only in
#' \code{PedCompare}).
#'   \item \code{FalsePos}(itives): no parent in reference pedigree, but
#' one was assigned based on the simulated data (P2only)
#'   \item \code{Mismatch}: different parents between the pedigrees
#'   }}
#' \item{cat}{Category of individual + parent, as a two-letter code where the
#'   first letter indicates the focal individual and the second the parent;
#'   G=Genotyped, D=Dummy, T=Total}
#' \item{parent}{dam or sire}
#'
#'
#' @seealso \code{\link{SimGeno}, \link{sequoia}, \link{PedCompare}}.
#'
#' @importFrom plyr adply
#'
#' @examples
#' # estimate proportion of parents that are genotyped (= 1 - ParMis)
#' prop_parents_genotyped <- c(
#'   dam = mean(unique(SeqOUT_griffin$Pedigree$dam) %in% rownames(Geno_griffin)),
#' sire = mean(unique(SeqOUT_griffin$Pedigree$sire) %in% rownames(Geno_griffin))
#' )
#'
#' # Example for parentage assignment only
#' conf_grif <- EstConf(Pedigree = SeqOUT_griffin$Pedigree,
#'                LifeHistData = SeqOUT_griffin$LifeHist,
#'                args.sim = list(nSnp = 150,   # no. in actual data, or what-if
#'                                SnpError = 5e-3,  # best estimate, or what-if
#'                                CallRate=0.9,     # from SnpStats()
#'                                ParMis=c(0.28, 0.22)),  # calc'd above
#'                args.seq = list(Err=5e-3, Module="par"),  # as in real run
#'                nSim = 1,   # try-out, proper run >=20 (10 if huge pedigree)
#'                nCores=1)
#'
#' # parent-pair confidence, per category (Genotyped/Dummy/None)
#' conf_grif$ConfProb
#'
#' # Proportion of true parents that was correctly assigned
#' 1 - apply(conf_grif$PedErrors, MARGIN=c('cat','parent'), FUN=sum, na.rm=TRUE)
#'
#' # add columns with confidence probabilities to pedigree
#' # first add columns with category (G/D/X)
#' Ped.withConf <- getAssignCat(Pedigree = SeqOUT_griffin$Pedigree,
#'                              SNPd = SeqOUT_griffin$PedigreePar$id)
#' Ped.withConf <- merge(Ped.withConf, conf_grif$ConfProb, all.x=TRUE,
#'                       sort=FALSE)  # (note: merge() messes up column order)
#' head(Ped.withConf[Ped.withConf$dam.cat=="G", ])
#'
#' # save output summary
#' \dontrun{
#' conf_griff[['Note']] <- 'You could add a note'
#' saveRDS(conf_grif[c('ConfProb','PedComp.fwd','RunParams','RunTime','Note')],
#'    file = 'conf_200SNPs_Err005_Callrate80.RDS')
#' }
#'
#' ## overall assignment rate (AR), error rate (ER) & runtime
#' AR_max <- sum(!is.na(Ped_griffin$dam)) + sum(!is.na(Ped_griffin$sire))
#' ER_max <- 2*nrow(Ped_griffin)
#' PCT <- conf_grif$PedComp.fwd[,'TT',,]   # Total-Total counts
#' list(AR = mean(apply(PCT[,'Match',],1,sum)/AR_max),  # sum over dam+sire
#'      ER = mean(apply(PCT[,c('Mismatch','P2only'),],1,sum)/ER_max),
#'      Time = mean(conf_grif$RunTime)/60)   # runtime in seconds --> minutes
#'
#'
#' ## P(actual FS | inferred as FS) etc.
#' \dontrun{
#' PairL <- list()
#' for (i in 1:length(conf_grif$Pedigree.inferred)) {  # nSim
#'   cat(i, "\t")
#'   PairL[[i]] <- ComparePairs(conf_grif$Pedigree.reference,
#'                              conf_grif$Pedigree.inferred[[i]],
#'                              GenBack=1, patmat=TRUE, ExcludeDummies = TRUE,
#'                              Return="Counts")
#' }
#' # P(actual relationship (Ped1) | inferred relationship (Ped2))
#' PairRel.prop.A <- plyr::laply(PairL, function(M)
#'                      sweep(M, MARGIN='Ped2', STATS=colSums(M), FUN="/"))
#' PairRel.prop <- apply(PairRel.prop.A, 2:3, mean, na.rm=TRUE) #avg across sims
#' round(PairRel.prop, 3)
#' # or: P(inferred relationship | actual relationship)
#' PairRel.prop2 <- plyr::laply(PairL, function(M)
#'    sweep(M, MARGIN='Ped1', STATS=rowSums(M), FUN="/"))
#' }
#'
#' \dontrun{
#' # confidence probability vs. sibship size
#' source('https://raw.githubusercontent.com/JiscaH/sequoiaExtra/main/conf_vs_sibsize.R')
#' conf_grif_nOff <- Conf_by_nOff(conf_grif)
#' conf_grif_nOff['conf',,'GD',]
#' conf_grif_nOff['N',,'GD',]
#' }
#'
#' @export


EstConf <- function(Pedigree = NULL,
                    LifeHistData = NULL,
                    args.sim = list(nSnp = 400, SnpError = 1e-3, ParMis=c(0.4, 0.4)),
                    args.seq = list(Module="ped", Err=1e-3, Tassign=0.5, CalcLLR = FALSE),
                    nSim = 10,
                    nCores = 1,
                    quiet=TRUE)
{

  # check input ----
  if (is.null(Pedigree))  stop("Please provide Pedigree")
  if (is.null(LifeHistData))  cli::cli_alert_warning("Running without `LifeHistData`")
  if (!is.null(args.sim) & !is.list(args.sim))  stop("args.sim should be a list or NULL")
  if (!is.null(args.seq) & !is.list(args.seq))  stop("args.seq should be a list or NULL")
  if (!is.wholenumber(nSim) || nSim<1 || length(nSim)>1)
    stop("nSim must be a single positive number")

  if (!(isTRUE(quiet) | isFALSE(quiet) | quiet=='very'))
    stop("'quiet' must be TRUE, FALSE, or 'very'")
  quiet.EC <- ifelse(quiet == "very", TRUE, FALSE)
  quiet <- ifelse(quiet %in% c("very", TRUE), TRUE, FALSE)
  if (!"quiet" %in% names(args.sim))  args.sim <- c(args.sim, list(quiet = quiet))
  if (!"quiet" %in% names(args.seq))  args.seq <- c(args.seq, list(quiet = quiet))

  if ("Err" %in% names(args.sim)) {
    args.sim[["SnpError"]] <- args.sim[["Err"]]
    args.sim[["Err"]] <- NULL    # common confusion, otherwise fuzy matching with 'ErrorFM'.
  }

  Ped.ref <- PedPolish(Pedigree, KeepAllColumns=FALSE)
  if (any(substr(unlist(Ped.ref),1,6) %in% c("sim_F0", "sim_M0"))) {
    stop("Please don't use 'sim_F' or 'sim_M' in reference pedigree")
  }

  if ("Module" %in% names(args.seq)) {
    ParSib <- ifelse(args.seq$Module == "ped", "sib", "par")
  } else {
    ParSib <- "sib"   # default Module = "ped"
  }

  if (!quiet.EC) {
    if (ParSib == "par") {
      cli::cli_alert_info("Simulating parentage assignment only ...")
    } else {
      cli::cli_alert_info("Simulating full pedigree reconstruction ...")
    }
  }


  # no. cores ----
  if (!is.null(nCores) && (!is.wholenumber(nCores) || nCores<1))
    stop("nCores must be a positive number, or NULL")

  if (is.null(nCores) || nCores>1) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      if (interactive() & !quiet.EC) {
        cli::cli_alert_info("Installing package {.pkg parallel} to run on multiple cores... ")
      }
      utils::install.packages("parallel")
    }
    maxCores <- parallel::detectCores()
    if (is.null(nCores)) {
      nCores <- maxCores -1
    } else if (nCores > maxCores) {
      nCores <- maxCores
      cli::cli_alert_warning("Reducing {.code nCores} to {maxCores}, as that's all you have")
    }
    if (nCores > nSim) {
      nCores <- nSim
    }
    if (!quiet.EC)  cli::cli_alert_info("Using {nCores} out of {maxCores} cores")
  }


  utils::flush.console()

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # function to simulate genotypes & infer pedigree ----

  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  SimInfer <- function(i, RefPedigree, args.sim,
                       LifeHistData, args.seq,
                       quiet.EC, ParSib)
  {
    if (!quiet.EC)  cat("i=", i, "\t", format(Sys.time(), "%H:%M:%S"), "\n")
    # nothing printed by parallel::parSapply (?)
    OUT.i <- list()
    GM <- do.call(sequoia::SimGeno, c(list(Pedigree=RefPedigree), args.sim))
    OUT.i$SimSNPd <- rownames(GM)
    OUT.i$RunTime <- system.time(Seq.i <- do.call(sequoia::sequoia,
                                                  c(list(GenoM = GM,
                                                         LifeHistData = LifeHistData,
                                                         DummyPrefix = c("sim_F", "sim_M"),
                                                         Plot = FALSE),
                                                    args.seq) ))["elapsed"]
    rm(GM)
    gc()
    if (ParSib == "par") {
      OUT.i$Pedigree.inferred <- Seq.i[["PedigreePar"]]
    } else {
      OUT.i$Pedigree.inferred <- Seq.i[["Pedigree"]]
    }
    return( OUT.i )
  }
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  # call above function, for 1 or >1 cores ----
  if (nCores>1) {
    cl <- parallel::makeCluster(nCores)
    AllOUT <- parallel::parLapply(cl, X=seq.int(nSim), fun=SimInfer,
                                  RefPedigree = Ped.ref,
                                  args.sim, LifeHistData, args.seq,
                                  quiet.EC=TRUE, ParSib)
    #                                  chunk.size=1)
    parallel::stopCluster(cl)
  } else {
    AllOUT <- plyr::llply(seq.int(nSim), .fun=SimInfer,
                          RefPedigree = Ped.ref,
                          args.sim, LifeHistData, args.seq,
                          quiet.EC, ParSib)
  }
  RunTime <- sapply(AllOUT, "[[", "RunTime")
  Pedigree.inferred <- sapply(AllOUT, "[", "Pedigree.inferred")
  SimSNPd <- sapply(AllOUT, "[", "SimSNPd")


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # confidence probabilities ----

  # counts for confidence probs, using PedCompare(Ped1 = Ped.inferred, Ped2 = Ped.ref).
  # function Count_conf() defined below in this file.
  CountsL <- list()
  for (i in 1:nSim) {
    CountsL[[i]] <- Count_conf(Pedigree.inferred[[i]], Ped.ref, SimSNPd[[i]])
  }

  # sum across iterations & calculate proportions
  # function Counts2Conf() defined below in this file.
  Conf.df <- Counts2Conf(CountsL)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # assignment errors (ignores co-parent) ----
  nSimz <- ifelse(nSim>1, nSim,2)  # else problems w R auto-dropping dimension
  PedComp.fwd <-  array(0, dim=c(nSimz, 7,5,2),
                        dimnames = list(iter = seq_len(nSimz),
                                        cat = c("GG", "GD", "GT", "DG", "DD", "DT", "TT"),
                                        class = c("Total", "Match", "Mismatch", "P1only", "P2only"),
                                        parent = c("dam", "sire")))
  for (i in 1:nSim) {
    PedComp.fwd[i,,,] <- PedCompare(Ped1 = Ped.ref,
                                    Ped2 = Pedigree.inferred[[i]],
                                    SNPd = SimSNPd[[i]],
                                    Symmetrical=FALSE, Plot=FALSE)$Counts
  }

  PedComp.tmp <- apply(PedComp.fwd, 2:4, sum)
  PedErrors <- sweep(PedComp.tmp[,c("P1only", "P2only","Mismatch"),], c(1,3),
                     PedComp.tmp[,"Total",], "/")
  PedErrors[c("GG", "GD", "DG", "DD"), "P2only", ] <- NA  # if parent in ref. pedigree, by def not P2only
  dimnames(PedErrors)[['class']] <- c("FalseNeg", "FalsePos", "Mismatch")



  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # out ----
  RunParams <- list(EstConf = namedlist(args.sim, args.seq, nSim, nCores),
                    SimGeno_default = formals(SimGeno),
                    sequoia_default = formals(sequoia),
                    sequoia_version = as.character(utils::packageVersion("sequoia")))

  return( list(ConfProb = Conf.df,
               PedErrors = PedErrors,
               Pedigree.reference = Ped.ref,
               LifeHistData = LifeHistData,
               Pedigree.inferred = Pedigree.inferred,
               SimSNPd = SimSNPd,
               PedComp.fwd = PedComp.fwd,
               RunParams = RunParams,
               RunTime = RunTime) )
}





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# counts for confidence calculations:
# Probability that an parent in the inferred pedigree is correct (is also the
# parent in the reference pedigree).  This function is called by EstConf(),
# which combines its result across multiple inferred pedigrees.

#' @keywords internal 
#' @noRd

Count_conf <- function(Ped.inferred, Ped.ref, SNPd) {

  # run PedCompare 'reversed': for each inferred category, you want the proportion correct
  PCCD <- PedCompare(Ped1 = Ped.inferred, Ped2 = Ped.ref, SNPd = SNPd,
                    Symmetrical=FALSE, Plot=FALSE)$Counts.detail

  OK <- list('G' = 'Match',
             'D' = 'Match',
             'X' = c('P2only', '_'))   # no parent assigned in Ped.inferred (includes non-assignable)

  conf.counts <- array(0, dim = c(4,3,3,3),
                 dimnames = c(list(c('N', paste0(c('dam', 'sire', 'pair'), '.OK'))),
                              dimnames(PCCD)[1:3]))

  # Totals per category (=denominator for the proportion)
  conf.counts['N',,,] <- apply(PCCD, c('id.cat', 'dam.cat', 'sire.cat'), sum)

  for (i in c('G','D','X')) {
    conf.counts['dam.OK' ,,i,] <- apply(PCCD[,i,,OK[[i]],], c('id.cat', 'sire.cat'), sum)
    conf.counts['sire.OK',,,i] <- apply(PCCD[,,i,,OK[[i]]], c('id.cat', 'dam.cat'), sum)
    for (j in c('G','D','X')) {
      conf.counts['pair.OK',,i,j] <- apply(PCCD[,i,j,OK[[i]],OK[[j]], drop=FALSE], 'id.cat', sum)
    }
  }

  return( conf.counts )
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# take list of outputs from Count_conf(), calculate proportions, output as dataframe

#' @keywords internal 
#' @noRd

Counts2Conf <- function(countsL) {

  # sum across iterations
  if (length(countsL)>1) {
    countsA <- apply( plyr::laply(countsL, function(x) x), 2:5, sum)
  } else {
    countsA <- countsL[[1]]
  }

  # array -> data.frame
  Conf.df <- plyr::adply(countsA[,c('G','D'),,], .margins=2:4)

  # counts -> proportions
  for (x in c('dam', 'sire', 'pair')) {
    Conf.df[[paste0(x,'.conf')]] <- Conf.df[[paste0(x,'.OK')]] / Conf.df[['N']]
    Conf.df[[paste0(x,'.OK')]] <- NULL   # drop column
  }

  # if parentage assignment only, drop all dummy related rows
  Ntot <- apply(countsA, 'id.cat', sum)
  if (Ntot[['D']] == 0) {
    Conf.df <- Conf.df[Conf.df$id.cat == 'G' & Conf.df$dam.cat %in% c('G','X') &
                         Conf.df$sire.cat %in% c('G','X'), ]
  }

  # if category='X', set respective conf to N/A
  no_dam <- Conf.df[['dam.cat']] == 'X'
  no_sire <- Conf.df[['sire.cat']] == 'X'
  Conf.df[['dam.conf']][no_dam] <- NA
  Conf.df[['sire.conf']][no_sire] <- NA
  Conf.df[['pair.conf']][no_dam | no_sire] <- NA

  # sort rows & columns
  Conf.df <- Conf.df[order(Conf.df$id.cat, Conf.df$dam.cat, Conf.df$sire.cat),
                     c('id.cat', 'dam.cat', 'sire.cat', 'dam.conf', 'sire.conf', 'pair.conf', 'N')]
  rownames(Conf.df) <- NULL

  return( Conf.df )
}
