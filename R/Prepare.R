
#=======================================================================
#=======================================================================
#' @title Specs to PARAM
#'
#' @description Convert 1-row dataframe \code{Specs} into list \code{PARAM},
#'   optionally including various other objects in the list.
#'
#' @param Specs 1-row dataframe, element of \code{\link{sequoia}} output list.
#' @param ... other objects to append to the list, such as \code{ErrM} and
#'   \code{quiet}.
#'
#' @return  A named list.
#'
#' @keywords internal

SpecsToParam <- function(Specs, ErrM = NULL, ErrFlavour = NULL,
                         dimGeno=NULL, ...)
{
  other_args <- list(...)
  if ("Module" %in% names(other_args)) {
    Specs$Module <- other_args$Module
  } else if (!"Module" %in% names(Specs)) {
    Specs$Module <- 'ped'  # default
  }
  if (!"Herm" %in% names(Specs)) {
    Specs$Herm <- switch(as.character(Specs$Complexity),
                         mono = "no",
                         simp = "no",
                         full = "no",
                         herm = "A",
                         herm2 = "B",
                         herm_simp = "A",
                         herm2_simp = "B")
    Specs$Complexity <- switch(as.character(Specs$Complexity),
                            mono = "mono",
                            simp = "simp",
                            full = "full",
                            herm = "full",
                            herm2 = "full",
                            herm_simp = "simp",
                            herm2_simp = "simp")
  }
  if (is.null(dimGeno)) {
    dimGeno <- with(Specs, c(NumberIndivGenotyped, NumberSnps))
  }  # else take from input (i.e., fresh GenoM)

  L <- with(Specs, namedlist(dimGeno,
                             Err = GenotypingErrorRate,
                             Tfilter,
                             Tassign,
                             nAgeClasses,
                             MaxSibshipSize,
                             Module,
                             DummyPrefix = c(DummyPrefixFemale,
                                             DummyPrefixMale),
                             Complex = Complexity,
                             Herm,
                             UseAge,
                             CalcLLR))

  L$ErrFlavour <- ifelse("ErrFlavour" %in% names(Specs),  # in $Specs from version 2.1 onwards
                         Specs$ErrFlavour,
                         ErrFlavour)
  if (!is.null(ErrM)) {
    L$ErrM <- ErrToM(ErrM, Return = "matrix")  # performs checks
  } else {
    L$ErrM <- ErrToM(L$Err, flavour = L$ErrFlavour, Return = "matrix")
  }

  if ("MaxMismatchOH" %in% names(Specs)) {  # DUP/OH/ME from version 2.0 onwards
    L$MaxMismatchV <- with(Specs, c(DUP = MaxMismatchDUP,
                                    OH = MaxMismatchOH,
                                    ME = MaxMismatchME))
  }

  if ("DummyPrefixHerm" %in% names(Specs))
    L$DummyPrefix <- c(L$DummyPrefix, Specs$DummyPrefixHerm)

  # note: code below is namedlist(), but not sure if/how to pass along names to other function with '...'
  L2 <- list(...)
  L2Names <- as.character(as.list( match.call())[-(1:5)])  # 1st is function name
  if (is.null(names(L2))) {  # all elements unnamed
    names(L2) <- L2Names
  } else {  # some elements unnamed
    names(L2) <- ifelse(names(L2)=="", L2Names, names(L2))
  }

  return( c(L, L2))
}


#=======================================================================
#=======================================================================

#' @title PARAM to Specs
#'
#' @description Convert list \code{PARAM} into 1-row dataframe \code{Specs}.
#'   Only to be called by \code{\link{sequoia}}.
#'
#' @param PARAM list with input parameters.
#' @param TimeStart  time at which \code{\link{sequoia}} run was started.
#' @param ErrFlavour  character name or function.
#'
#' @return  The 1-row \code{Specs} dataframe.
#'
#' @keywords internal

ParamToSpecs <- function(PARAM, TimeStart, ErrFlavour)
{

  DF <- with(PARAM,
             data.frame(NumberIndivGenotyped = dimGeno[1],
                        NumberSnps = dimGeno[2],
                        GenotypingErrorRate = ifelse(length(Err)==1,
                                                     Err,
                                                     signif(1 - mean(c(diag(ErrM),
                                                                       ErrM[2,2])), 3)),  # under HWE, MAF=0.5: 0.25-0.5-0.25
                        MaxMismatchDUP = MaxMismatchV["DUP"],
                        MaxMismatchOH = MaxMismatchV["OH"],
                        MaxMismatchME = MaxMismatchV["ME"],
                        Tfilter = Tfilter,
                        Tassign = Tassign,
                        nAgeClasses = nAgeClasses,
                        MaxSibshipSize = MaxSibshipSize,
                        Module = as.character(Module),
                        DummyPrefixFemale = DummyPrefix[1],
                        DummyPrefixMale = DummyPrefix[2],
                        DummyPrefixHerm = DummyPrefix[3],
                        Complexity = Complex,
                        Herm = Herm,
                        UseAge = UseAge,
                        CalcLLR = CalcLLR,
                        ErrFlavour = ifelse(length(Err)==9,
                                            "customMatrix",
                                            ifelse(is.function(ErrFlavour),
                                                   "customFunction",
                                                   as.character(ErrFlavour))),
                        SequoiaVersion = as.character(utils::packageVersion("sequoia")),
                        TimeStart = TimeStart,
                        TimeEnd = Sys.time(),   # this function is called at the end of sequoia()
                        stringsAsFactors = FALSE))
  if (PARAM$Herm == "no")  DF <- DF[, colnames(DF)!="DummyPrefixHerm"]
  rownames(DF) <- "Specs"
  return( DF )
}


#=======================================================================
#=======================================================================
#' @title PARAM to FortPARAM
#'
#' @description Convert list \code{PARAM} into a list with integer-only and
#'   double-only vectors, to be passed to Fortran.
#'
#' @param PARAM list with input parameters.
#' @param fun  function from which \code{\link{MkFortParams}} is called,
#'   determines which elements are included in the output list.
#'
#' @return  A list with elements
#'   \item{Ng}{Integer, number of individuals}
#'   \item{SpecsInt}{8 integers:
#'     \itemize{
#'       \item nSnp
#'       \item MaxMismatchV; DUP - OH - ME
#'       \item MaxSibshipSize
#'       \item Complx, 0=mono, 1=simp, 2=full
#'       \item quiet, -1=verbose, 0=FALSE, 1=TRUE
#'       \item nAgeCl, nrow(AgePriors)
#'     }
#'   }
#'   \item{SpecsDbl}{2 double precision numbers:
#'     \itemize{
#'       \item Tfilter  (< 0)
#'       \item Tassign  (> 0)
#'     }
#'   }
#'   \item{ErrM}{double, 3x3 matrix passed as length-9 vector}
#'   \item{SpecsIntMkPed}{\code{fun='main'} only
#'     \itemize{
#'       \item AgeEffect, 0=no, 1=yes, 2=extra
#'       \item CalcLLR, 0=FALSE, 1=TRUE
#'       \item Herm, 0=no, 1= dam/sire distinction, 2=no dam/sire distinction
#'     }
#'   }
#'   \item{SpecsIntAmb}{\code{fun='mayberel'} only
#'     \itemize{
#'       \item ParSib  1=par, 2=ped
#'       \item nAmbMax
#'     }
#'   }
#'
#' @keywords internal

MkFortParams <- function(PARAM, fun="main")
{
  FP <- list()
  FP$Ng = as.integer(PARAM$dimGeno[1])
  FP$SpecsInt <- with(PARAM,
                      c(nSnp = dimGeno[2],      # 1
                        MaxMismatchV,           # 2 - 4
                        MaxSibshipSize = ifelse(exists("MaxSibshipSize"),
                                                MaxSibshipSize,
                                                100),  # 5
                        Complx = switch(as.character(Complex),       # 6
                                        mono = 0,
                                        simp = 1,
                                        full = 2),
                        quiet = ifelse(quiet == "verbose", -1,     # 7 FALSE=0, TRUE=1
                                       ifelse(quiet == "very", 1,
                                              as.integer(as.logical(quiet)))),
                        nAgeCl = nAgeClasses,      # 8
                        Herm = switch(as.character(Herm),   # 9
                                      no = 0,
                                      A = 1,
                                      B = 2) ))
  FP$SpecsDbl <- as.double(c(PARAM$Tfilter,
                             PARAM$Tassign))
  FP$ErrM <- as.double(PARAM$ErrM)

  if (fun == "main" | fun == "CalcOH") {
    # ParSib added in SeqParSib()
    FP$SpecsIntMkPed <- with(PARAM, c(AgeEffect = ifelse(exists("UseAge"),
                                                     switch(UseAge,
                                                         no = 0,
                                                         yes = 1,
                                                         extra = 2),
                                                     1),
                                      CalcLLR = as.logical(CalcLLR) ))
  } else if (fun == "mayberel") {
    FP$SpecsIntAmb <- with(PARAM, c(ParSib = switch(Module,
                                                    par = 1,
                                                    ped = 2),
                                    nAmbMax = MaxPairs))
    #  } else if (fun == "CalcOH") {
    #    # do nothing
  } else if (fun ==  "CalcPairs") {
    # do nothing
  } else {
    stop("invalid 'fun'")
  }

  return( FP )
}


#=======================================================================
#=======================================================================
#' @title Check if input parameters are valid
#'
#' @description Check if input parameter value is of the proper kind and value.
#'
#' @param PARAM list with input parameters
#'
#' @return  Nothing except errors, warnings, and messages
#'
#' @keywords internal

CheckParams <- function(PARAM)
{

  # checking function ----
  chk <- function(x, type, valid, n=1, allowNeg = FALSE) {
    if (!x %in% names(PARAM))  return()
    if (!type %in% c("value", "int", "dbl"))  stop("type not supported")

    ErrMsg <- switch(type,
                     value = paste(sQuote(valid), collapse=" or "),
                     int = paste0(ifelse(n==1, "a", n),
                                  ifelse(allowNeg, " ", " positive "),
                                  "whole number", ifelse(n==1, "", "s")),
                     dbl = paste(ifelse(allowNeg, "", " positive"), "number"))
    ErrMsg <- paste(sQuote(x), "must be", ErrMsg)

    xx <- PARAM[[x]]
    if (is.null(xx) || length(xx) != n) {
      stop(ErrMsg, call.=FALSE)
    } else if (type == "value" && !xx %in% valid) {
      stop(ErrMsg, call.=FALSE)
    } else if (type == "int" && (any(!is.wholenumber(xx)) ||
                                 (!allowNeg & any(xx < 0, na.rm=TRUE)))) {
      stop(ErrMsg, call.=FALSE)
    } else if (type=="dbl" && (any(!is.double(xx)) ||
                               (!allowNeg & any(xx < 0, na.rm=TRUE)))){
      stop(ErrMsg, call.=FALSE)
    }
  }  # TODO: use match.arg() ; stopifnot(). See https://blog.r-hub.io/2022/03/10/input-checking/

  # check if in PARAM list ----
  chk("quiet", "value", valid = c(TRUE, FALSE, 'verbose', 'very'))
  chk("CalcLLR", "value", valid = c(TRUE, FALSE))
  chk("GenoProb", "value", valid = c(TRUE, FALSE))
  chk("Plot",  "value", valid = c(TRUE, FALSE))
  chk("dropPar", "value", valid = c(0, 1, 2, 12))
  chk("Module", "value", valid=c("pre", "dup", "par", "ped"))
  chk("ParSib", "value", valid=c("dup", "par", "sib"))
  chk("Complex", "value", valid=c("full", "simp", "mono"))  # 'herm'
  chk("Herm", "value", valid=c("no", "A", "B"))
  chk("UseAge", "value", valid=c("yes", "no", "extra"))

  chk("MaxMismatchV", "int", n=3)
  chk("MaxSibshipSize", "int")
  chk("nAgeClasses", "int")
  chk("MaxPairs", "int")

  chk("Tfilter", "dbl", allowNeg = TRUE)
  chk("Tassign", "dbl")


  # other
  if ("DummyPrefix" %in% names(PARAM)) {
    if (!length(PARAM$DummyPrefix) %in% c(2,3) |
        any(make.names(PARAM$DummyPrefix) != PARAM$DummyPrefix))
      stop("`DummyPrefix` should be a length-2 character vector with syntactically valid names",
           call.=FALSE)
  }
  if ("MaxMismatch" %in% names(PARAM)) {
    if (!is.null(PARAM$MaxMismatch) && !is.na(PARAM$MaxMismatch))
      cli::cli_alert_info("NOTE: `MaxMismatch` is deprecated & ignored;",
              "now calculated internally by `CalcMaxMismatch()`")
  }
}


#=======================================================================
#=======================================================================
#' @title check AgePrior
#'
#' @description Check that the provided AgePrior is in the correct format
#'
#' @param AgePrior matrix with 'MaxAgeParent' rows and columns M-P-FS-MHS-PHS
#'
#' @return AgePrior in corrected format, if necessary
#'
#' @keywords internal

CheckAP <- function(AgePrior) {
  if (is.null(AgePrior)) {
    stop("Please provide AgePrior")
  } else if (is.data.frame(AgePrior)) {
    AgePrior <- as.matrix(AgePrior)
  } else if (!is.matrix(AgePrior)) {
    stop("AgePrior must be a matrix")
  }

  if (any(AgePrior < 0 | AgePrior > 1000) | any(!is.double(AgePrior)))
    stop("AP must be a numeric matrix with values >0 & <1000")

  if (!all(c("M", "P", "FS", "MS", "PS") %in% colnames(AgePrior)))
    stop("AgePrior must have columns M, P, FS, MS, PS")
  AgePrior <- AgePrior[, c("M", "P", "FS", "MS", "PS")]

  if (any(as.numeric(rownames(AgePrior)) < 0))
    AgePrior <- AgePrior[as.numeric(rownames(AgePrior)) >= 0, ]

  if (!all(apply(AgePrior, 2, function(x) any(x > 0))))
    stop("AgePriors error: some relationships are impossible for all age differences")

  span <- apply(AgePrior[,c("M","P")], 2, function(A) diff(range(which(A>0))))
  span["F"] <- diff(range(which(AgePrior[,"M"]>0 & AgePrior[,"P"]>0)))
  for (r in c("M", "P", "F")) {
    if(any(any(AgePrior[(span[r]+2):nrow(AgePrior), paste0(r,"S")] > 0)))  # 1st row = agedif 0
      stop(paste("Sibling column ", paste0(r,"S"), " not consistent with parent column"))
  }

  return( AgePrior )
}



#=======================================================================
#=======================================================================

#' @title Check and recode mtSame matrix
#'
#' @description Recode to 1=different, 0=same
#'
#' @param mtSame matrix indicating whether individuals (might) have the same
#'   mitochondrial haplotype (1),  or definitely not (0). Not all individuals
#'   need to be included and order is not important, may not be square.
#' @param gID rownames of `GenoM`
#'
#' @return  square gID x gID matrix indicating whether individuals have
#'   definitely a different mitochondrial haplotype (1), or (possibly) the same
#'   (0).
#'
#' @keywords internal

mtSame2Dif <- function(mtSame = NULL, gID = NULL)
{
  if (is.null(gID))  stop('must provide gID')

  if (is.null(mtSame)) {
    mtDif <- matrix(0, length(gID), length(gID), dimnames=list(gID,gID))
    return( mtDif )
  }

  if (!is.matrix(mtSame))  stop('mtSame must be a 0/1 matrix or NULL', call.=FALSE)
  if (!all(mtSame %in% c(0,1,NA)))   stop('mtSame must be a 0/1 matrix or NULL', call.=FALSE)
  if (!any(rownames(mtSame) %in% gID) | !any(colnames(mtSame) %in% gID)) {
    stop('dimnames of mtSame do not match IDs in GenoM', call.=FALSE)
  }

  mtDif <- 1 - inflate(mtSame[rownames(mtSame) %in% gID, colnames(mtSame) %in% gID], gID, na=1)
  mtDif[is.na(mtDif)] <- 0
  return( mtDif )
}
