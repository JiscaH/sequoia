#' @title Convert Genotype Data
#'
#' @description Convert genotype data in various formats to sequoia's
#'   1-column-per-marker format, PLINK's ped format, or Colony's
#'   2-columns-per-marker format.
#'
#' @param InFile character string with name of genotype file to be converted.
#' @param InData dataframe, matrix or
#'   \code{\link[adegenet:genlight-class]{genlight}} object with genotypes to be
#'   converted.
#' @param InFormat One of 'seq' (sequoia), 'ped' (PLINK .ped file), 'col'
#'   (COLONY), 'raw' (PLINK --recodeA), 'vcf' (requires library \code{{vcfR}}),
#'   'single' (1 column per SNP), or 'double' (2 columns per SNP); see Details.
#' @param OutFile character string with name of converted file. If NA, return
#'   matrix with genotypes in console (default); if NULL, write to
#'   'GenoForSequoia.txt' in current working directory.
#' @param OutFormat as \code{InFormat}; only 'seq', 'col', and 'ped' are
#'   implemented. For 'ped' also a sham .map file is created, so that the file
#'   can be read by PLINK. Only for 'ped' are extensions .ped & .map added to
#'   the specified OutFile filename.
#' @param Missing vector with symbols interpreted as missing data. '0' is
#'   missing data for \code{InFormat}s 'col' and 'ped' only.
#' @param sep vector with field separator strings that will be tried on
#'   \code{InFile}. Ignored if package \pkg{data.table} is present or if
#'   \code{InFormat='vcf'} or 'vcf.gz'. The \code{OutFile} separator uses the
#'   \code{\link[utils]{write.table}} default, i.e. one blank space.
#' @param header a logical value indicating whether the file contains a header
#'   as its first line. If NA (default), set to TRUE for 'raw', and FALSE
#'   otherwise.
#' @param IDcol number giving the column with individual IDs; 0 indicates the
#'   rownames (for InData only). If NA (default), set to 2 for InFormat 'raw'
#'   and 'ped', and otherwise to 1 for InFile and 0 (rownames) for InData,
#'   except when InData has a column labeled 'ID'.
#' @param FIDcol  column with the family IDs, if any are wished to
#'   be used. This is column 1 for InFormat 'raw' and 'seq', but those are by
#'   default not used.
#' @param FIDsep string used to paste FID and IID together into a composite-ID
#'   (value passed to \code{\link{paste}}'s \code{collapse}). This joining can
#'   be reversed using \code{\link{PedStripFID}}.
#' @param dropcol  columns to exclude from the output data, on top of IDcol and
#'   FIDcol (which become rownames). When NA, defaults to columns 3-6 for
#'   InFormat 'raw' and 'seq'. Can also be used to drop some SNPs, see example
#'   below on how to do this for the 2-columns-per-SNP input formats.
#' @param quiet suppress messages and warnings.
#'
#' @return A genotype matrix in the specified output format; the default sequoia
#'   format ('seq') has 1 column per SNP coded in 0/1/2 format (major homozygote
#'   /heterozygote /minor homozygote) with -9 for missing values, sample IDs in
#'   row names and SNP names in column names. If 'OutFile' is specified, the
#'   matrix is written to this file and nothing is returned inside R.
#'
#' @details The first two arguments are interchangeable, and can be given
#'   unnamed. The first argument is assumed to be a file name if it is of class
#'   'character' and length 1, and to be the genetic data if it is a matrix or
#'   dataframe.
#'
#'   If package \pkg{data.table} is detected, \code{\link[data.table]{fread}}
#'   is used to read in the data from file. Otherwise, a combination of
#'   \code{\link{readLines}} and \code{\link{strsplit}} is used.
#'
#'
#' @section Input formats:
#' The following formats can be specified by \code{InFormat}:
#' \describe{
#'   \item{seq}{(sequoia) genotypes are coded as 0, 1, 2, missing as \eqn{-9} (in input
#'  any negative number or NA are OK),
#'   in 1 column per marker. Column 1 contains IDs, there is no header row.}
#'   \item{ped}{(PLINK) genotypes are coded as A, C, T, G, missing as 0, in 2
#'   columns per marker. The first 6 columns are descriptive (1:FID, 2:IID, 3 to
#'   6 ignored). If an associated .map file exists, SNP names will be read from
#'   there.}
#'   \item{raw}{(PLINK) genotypes are coded as 0, 1, 2, missing as NA, in 1
#'   column per marker. The first 6 columns are descriptive (1:FID, 2:IID, 3 to
#'   6 ignored), and there is a header row. This is produced by PLINK's option
#'   --recodeA}
#'   \item{col}{(Colony) genotypes are coded as numeric values, missing as 0, in
#'   2 columns per marker. Column 1 contains IDs.}
#'   \item{vcf}{(VCF) genotypes are coded as '0/0','0/1','1/1', variable number
#'     of header rows followed by 1 row per SNP, with various columns of
#'     metadata followed by 1 column per individual. Requires package
#'     \pkg{vcfR}.}
#'   \item{single}{1 column per marker, otherwise unspecified}
#'   \item{double}{2 columns per marker, otherwise unspecified}
#'   }
#'  For each \code{InFormat}, its default values for \code{Missing, header,
#'  IDcol, FIDcol}, and \code{dropcol} can be overruled by specifying the
#'  corresponding input parameters.
#'
#' @section Error messages:
#'   Occasionally when reading in a file \code{GenoConvert} may give an error
#'   that 'rows have unequal length'. \code{GenoConvert} makes use of
#'   \code{\link{readLines}} and \code{\link{strsplit}}, which is much faster
#'   than \code{\link{read.table}} for large datafiles, but also more sensitive
#'   to unusual line endings, unusual end-of-file characters, or invisible
#'   characters (spaces or tabs) after the end of some lines. In these cases,
#'   try to read the data from file using read.table or read.csv, and then use
#'   \code{GenoConvert} on this dataframe or matrix, see example.
#'
#'  Any warnings generated by \code{\link{CheckGeno}} regarding SNPs scored for
#'  few individuals and/or individuals scored for few SNPs etc. are only for
#'  your information; none are excluded from \code{GenoConvert}'s output, but
#'  these SNPs and/or individuals will be excluded during pre-processing of the
#'  data in any of the other functions in this package.
#'
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{CheckGeno}, \link{SnpStats}, \link{LHConvert}}.
#'
#' @examples
#' \dontrun{
#' # Requires PLINK installed & in system PATH:
#'
#' # tinker with window size, window overlap and VIF to get a set of
#' # 400 - 800 markers (100-200 enough for just parentage):
#' system("cmd", input = "plink --file mydata --indep 50 5 2")
#' system("cmd", input = "plink --file mydata --extract plink.prune.in
#'   --recodeA --out PlinkOUT")
#'
#' GenoM <- GenoConvert(InFile = "PlinkOUT.raw", InFormat='raw')
#' # which is the same as
#' GenoM <- GenoConvert(PlinkOUT.raw, InFormat='single',
#'                     IDcol=2, dropcol=c(1,3:6), header=TRUE)
#' # (but it will complain that InFormat='single' is not consistent with .raw
#' # file extension)
#'
#' # save time on file conversion next time:
#' write.table(GenoM, file="Geno_sequoia.txt", quote=FALSE, col.names=FALSE)
#' GenoM <- as.matrix(read.table("Geno_sequoia.txt", row.names=1, header=FALSE))
#'
#' # drop some SNPs, e.g. after a warning of >2 alleles:
#' dropSNP <- c(5,68,101,128)
#' GenoM <- GenoConvert(ColonyFile, InFormat = "col",
#'                      dropcol = 1 + c(2*dropSNP-1, 2*dropSNP) )
#'
#' # circumvent a 'rows have unequal length' error:
#' GenoTmp <- as.matrix(read.table("mydata.txt", header=TRUE, row.names=1))
#' GenoM <- GenoConvert(InData=GenoTmp, InFormat="single", IDcol=0)
#'
#' # can also write to file, e.g. simulated genotypes:
#' GenoConvert(Geno_A, InFormat='seq', OutFormat='ped', OutFile = sim_genotypes)
#' }
#'
#' @importFrom stats na.exclude
#'
#' @export


GenoConvert <- function(InData = NULL,
                        InFile = NULL,
                        InFormat = "raw",
                        OutFile = NA,
                        OutFormat = "seq",
                        Missing = c("-9", "NA", "??", "?","NULL","-1",
                                    c("0")[InFormat %in% c("col","ped")]),
                        sep = c(" ", "\t", ",", ";"),
                        header = NA,
                        IDcol = NA, #
                        FIDcol = NA,
                        FIDsep = "__",
                        dropcol = NA,
                        quiet = FALSE) {

  if (!(isTRUE(quiet) | isFALSE(quiet)))  stop("'quiet' must be TRUE or FALSE")
  if (is.null(InFile) & is.null(InData)) {
    stop("please provide 'InFile' or 'InData'")

  } else if (!is.null(InFile) & !is.null(InData)) {
    stop("please provide either 'InFile' or 'InData', not both")
  }
  if (length(InData)==1 & inherits(InData, "character")) {
    InFile <- InData
    InData <- NULL
  } else if (is.matrix(InFile) | is.data.frame(InFile)) {
    InData <- InFile
    InFile <- NULL
  }

  if (!is.null(InFile)) {
    if (is.character(InFile)) {
      if (!file.exists(InFile)) stop("InFile ", InFile, ' not found')
    } else {
      stop("'InFile' in unknown format, should be character string.")
    }
  }

  if (!InFormat %in% c("raw", "seq", "col", "ped", "single", "double", "vcf")) {
    stop("invalid InFormat")
  }
  InExt <- tools::file_ext(InFile)
  if (InExt=='gz') {  # allow for .vcf.gz
    tmp <- tools::file_ext(gsub('.gz','', InFile))
    if (tmp=='vcf')  InExt <- 'vcf'
  }
  if (!is.null(InFile)) {
    for (xx in c('ped', 'raw', 'vcf', 'vcf.gz')) {
      if (InExt == xx & InFormat != xx) {
        if (InFormat == 'raw') {
          InFormat = xx   # change from default
        } else {
          stop("InFormat ", InFormat, " not consistent with file extension .",xx," for InFile")
        }
      } else if (InExt != xx & InFormat == xx) {
        stop("File extension ", InExt, " not consistent with InFormat ", xx)
      }
    }
  }
  if (OutFormat %in% c("raw", "single")) {
    stop("OutFormat not (yet) implemented")
  } else if (!OutFormat %in% c("seq", "col", "ped")) {
    stop("invalid OutFormat")
  }

  if (!quiet)  cli::cli_alert_info('Input format is: {InFormat}')

  if (is.na(header))
    header <- ifelse(InFormat == "raw", TRUE, FALSE)

  if (!is.na(FIDcol) & FIDsep %in% c("", " ", "\t", "\n")) stop("sep can not be whitespace")

  if (is.na(IDcol)) {
    IDcol <- ifelse("ID" %in% colnames(InData),
                    which(colnames(InData) == "ID"),
                    switch(InFormat,
                           ped = 2,
                           raw = 2,
                           vcf = 0,
                           ifelse(!is.null(InFile), 1, 0)))  # 1st column; rownames
  }
  if (InFormat %in% c("raw", "ped") & all(is.na(dropcol))) {
    dropcol <- c(3:6)
    if (is.na(FIDcol))  dropcol <- c(1, dropcol)
  }

  if (OutFormat == "seq" & is.null(OutFile)) {  # NA=to console, NULL=to file
    OutFile <- "GenoForSequoia.txt"

  } else if (is.null(OutFile)) {
    stop("please provide 'OutFile'")
  }
  if (interactive() & !quiet & !(is.na(OutFile))) {
    if (file.exists(OutFile)) {
      ANS <- readline(prompt = paste("WARNING: ", OutFile, " will be overwritten.",
                                     "\n Press <N> to abort, or any other key to continue."))
    } else {
      ANS <- readline(prompt = paste("Genotypes will be written to ", OutFile,
                                     " \n Press <N> to abort, or any other key to continue."))
    }
    if (substr(ANS, 1, 1) %in% c("N", "n")) stop()
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~
  SNPnames <- NULL

  if (!is.null(InData)) {
    GenoTmp <- as.matrix(InData)
    rm(InData)

  } else if (InFormat == 'vcf') {
    if (!requireNamespace("vcfR", quietly = TRUE)) {
      if (interactive() & !quiet) {
        ANS <- readline(prompt = cli::cli_text("package {.pkg vcfR} not found. Install Y/N? "))
        if (!substr(ANS, 1, 1) %in% c("Y", "y")) stop(call.=FALSE)
      }
      utils::install.packages("vcfR")
    }
    GenoTmp <- vcfR::read.vcfR(InFile)
    GenoTmp <- vcfR::vcfR2genlight(GenoTmp)
    GenoTmp <- as.matrix(GenoTmp)
    if (all(substr(rownames(GenoTmp),1,2) == '0_') & is.na(FIDcol)) {  # FID prefix added by PLINK
      rownames(GenoTmp) <- gsub('^0_', '', rownames(GenoTmp))
    }

  } else if (!is.null(InFile)) {
    if (requireNamespace("data.table", quietly = TRUE)) {
      GenoTmp <- data.table::fread(InFile, header=header)
      GenoTmp <- as.data.frame(GenoTmp)
    } else {
      GenoTmp <- readLines(InFile, warn=FALSE)
      if (header) {
        SNPnames <- GenoTmp[1]
        GenoTmp <- GenoTmp[-1]
      }
      # tidy up blanks  (NOTE: trimws() is very slow)
      GenoTmp <- lapply(GenoTmp, function(e) {
        e <- gsub("[[:blank:]]+", " ", e)
        e <- gsub('^ ','',e)
        e <- gsub(' $','',e) })
      GenoTmp <- sapply(GenoTmp, strsplit, split=sep[1])
      # check if separator worked; if not, try out different separators
      if (length(GenoTmp[[1]])==1) {
        for (s in sep[-1]) {
          GenoTmp <- strsplit(GenoTmp, split = s)
          if (length(GenoTmp[[1]]) > 1) break
        }
      }
      if (length(GenoTmp[[1]])==1) {
        stop("unknown column separator, expecting ' ' (space), '\t' (tab), ',' or ';'")
      }
      if (length(table(sapply(GenoTmp, length))) > 1) {
        stop("rows have unequal length")
      }
      GenoTmp <- plyr::ldply(GenoTmp)
    }
  }
  if (nrow(GenoTmp)<2)  stop("Genotype matrix must have at least 2 individuals")
  if (ncol(GenoTmp)<2)  stop("Genotype matrix must have at least 2 SNPs")

  #~~~~~~~~~
  # IDs
  if (IDcol==0) {
    IDs_geno <- rownames(GenoTmp)
  } else {
    IDs_geno <- trimws( GenoTmp[, IDcol] )
  }
  if (!is.na(FIDcol) & InFormat!='vcf') {
    FID <- GenoTmp[, FIDcol]
    IDs_geno <- paste(FID, IDs_geno, sep=FIDsep)
  }
  if (any(duplicated(IDs_geno))) {
    dup_msg <- paste("`GenoM` has duplicate IDs in ",
                     ifelse(IDcol==0, 'rownames', paste0('column ', IDcol)),": ")
    if (length(unique(IDs_geno))/length(IDs_geno) > 0.5) {
      dup_IDs <- utils::head(IDs_geno[duplicated(IDs_geno)],3)
      cli::cli_alert_danger(dup_msg,
                            cli::cli_li(c(dup_IDs, "...")))
    } else {
      cli::cli_alert_danger(dup_msg,
                            cli::cli_li(c(IDs_geno[1:2], "...", utils::tail(IDs_geno,1))))
    }
    stop("Please exclude or rename these samples, or specify `IDcol` or `FIDcol` (or set `InFormat`).")
  }




  #~~~~~~~~~
  # drop non-genotype columns
  dropcol <- na.exclude(c(FIDcol, IDcol, dropcol))
  if (any(dropcol!=0)) {
    if (!quiet)   cli::cli_alert_info('Dropping columns: {setdiff(dropcol,0)} from input')
    GenoTmp <- GenoTmp[, -dropcol]
  }
  GenoTmp <- as.matrix(GenoTmp)


  #~~~~~~~~~
  # SNP names
  if (!is.null(InFile) & InFormat=='ped') {
    # check if map file exists & read in SNP names from there
    MapFile <- gsub('.ped$', '.map', InFile)
    if (file.exists(MapFile))
      SNPnames <- utils::read.table(MapFile, header=FALSE)[,2]
  }
  if (is.null(SNPnames)) {
    nSnp <- ifelse(InFormat %in% c("ped", "col", "double"), ncol(GenoTmp)/2, ncol(GenoTmp))
    if (is.null(colnames(GenoTmp)) |
        all(colnames(GenoTmp) == paste0('V', seq_len(ncol(GenoTmp)))) |
        ncol(GenoTmp)!=nSnp) {
      SNPnames <- paste0('SNP', formatC(seq_len(nSnp), width=5, flag=0))
    } else {
      SNPnames <- colnames(GenoTmp)
    }
  }
  # remove suffixes from SNP names
  if (InFormat=='raw') {
    for (x in c('A','C','T','G'))  SNPnames <- gsub(paste0('_',x,'$'), '', SNPnames)
  }

  #~~~~~~~~~
  # Recode all missing values to NA
  for (misX in Missing) {
    if (any(GenoTmp == misX, na.rm=TRUE))
      GenoTmp <- array(gsub(misX, NA, GenoTmp, fixed=TRUE), dim=dim(GenoTmp))
  }

  if (InFormat %in% c("raw", "seq")) {
    if (!all(GenoTmp %in% c(0,1,2, NA)))
      stop(paste("Unexpected value! When InFormat=", InFormat, ", genotypes should be coded as 0/1/2.\n",
                 "Choose InFormat='single' and/or different Missing, or check data."))
  } else if (InFormat =="single") {
    if (!all(nchar(GenoTmp)==2, na.rm=T)) {
      if (all(GenoTmp %in% c(0,1,2, NA))) {
        InFormat <- 'seq'
      } else {
        stop(paste("Unexpected value! When InFormat='single', genotypes should be coded as 2 digits or characters.\n",
                   "Choose other InFormat, or check data."))
      }
    }
  }

  if (InFormat %in% c("col", "ped", "single", "double")) {  # A/C/T/G -> 0/1/2

    if (InFormat %in% c("ped", "col", "double")) {
      GCA <- array(dim=c(2, nrow(GenoTmp), ncol(GenoTmp)/2))
      GCA[1,,] <- GenoTmp[, seq(1,ncol(GenoTmp)-1,2)]
      GCA[2,,] <- GenoTmp[, seq(2,ncol(GenoTmp),2)]
    } else if (InFormat %in% c("single")){
      GCA <- array(dim=c(2, nrow(GenoTmp), ncol(GenoTmp)))
      GCA[1,,] <- substr(GenoTmp,1,1)
      GCA[2,,] <- substr(GenoTmp,2,2)
    }
    UniqueAlleles <- sort(na.exclude(unique(c(GCA))))
    Alleles <- apply(GCA, 3, function(x) table(factor(x, levels = UniqueAlleles)))

    if (is.matrix(Alleles)) {
      NumAlleles <- apply(Alleles, 2, function(x) sum(x>0))
      minorAllele <- apply(Alleles, 2, function(x) names(sort(x[x>0]))[1])
    } else {  # list
      NumAlleles <- sapply(Alleles, length)
      minorAllele <- sapply(Alleles, function(x)  names(sort(x))[1])
    }
    if (any(NumAlleles > 2)) {
      n.problem <- sum(NumAlleles > 2)
      stop(paste("There are", n.problem, "SNPs with >2 alleles  ",
                 ifelse(n.problem<=10, paste(which(sapply(Alleles, length) > 2), collapse="-"),"")))
    }
    if (any(NumAlleles ==1) & !quiet & OutFormat!='seq')
      cli::cli_alert_warning(paste("There are {sum((NumAlleles ==1))} monomorphic SNPs"))
    GenoTmp2 <- sapply(1:length(minorAllele), function(i) {
      apply(GCA[,,i], 2, function(x) sum(x == minorAllele[i]))  })

  } else {
    AllHom0 <- apply(GenoTmp, 2, function(x) all(na.exclude(x) == 0))
    AllHom2 <- apply(GenoTmp, 2, function(x) all(na.exclude(x) == 2))
    if ((any(AllHom0) | any(AllHom2)) & !quiet & OutFormat!='seq') {
      cli::cli_alert_warning(paste("There are {sum(AllHom0)+sum(AllHom2)} monomorphic SNPs"))
    }

    GenoTmp2 <- matrix(as.numeric(GenoTmp), nrow(GenoTmp))
  }

  rownames(GenoTmp2) <- IDs_geno
  colnames(GenoTmp2) <- SNPnames
  rm(GenoTmp)

  #~~~~~~~~~
  # Format output
  if (OutFormat == "seq") {
    GenoOUT <- GenoTmp2
    GenoOUT[is.na(GenoOUT)] <- -9
    CheckGeno(GenoOUT, quiet=quiet, Return = "excl")  # returns invisibly

  } else if (OutFormat %in% c("col", "ped")) {
    dc <- list("0" = c(1,1), "1" = c(1,2), "2" = c(2,2), "-9" = c(0,0))
    GenoTmp2[is.na(GenoTmp2)] <- -9
    GenoA <- array(dim=c(nrow(GenoTmp2), 2, ncol(GenoTmp2)))
    for (i in 1:nrow(GenoTmp2)) {
      GenoA[i,,] <- sapply(GenoTmp2[i,], function(z) dc[[as.character(z)]])
    }
    GenoOUT <- matrix(GenoA, nrow(GenoTmp2))

    if (OutFormat == "col") {
      row.names(GenoOUT) <- IDs_geno

    } else if (OutFormat == 'ped') {
      FamOUT <- data.frame(FID = 0,
                           IID = IDs_geno,
                           sire = 0,
                           dam = 0,
                           sex = 0,
                           pheno = 0)

      MapOUT <- data.frame(chrom = 1,
                           SNP = SNPnames,
                           pos_cM = 0,
                           pos_bp = 0)
    }
    #  } else {
    #    stop("OutFormat not implemented")  # caught above
  }

  if (!is.na(OutFile)) {  # output to file
    if (OutFormat == 'ped') {
      utils::write.table(cbind(FamOUT, GenoOUT),
                         file = paste0(OutFile, ".ped"),
                         row.names = FALSE, col.names = FALSE, quote = FALSE)
      utils::write.table(MapOUT, file = paste0(OutFile, ".map"),
                         row.names = FALSE, col.names = FALSE, quote = FALSE)

    } else {
      utils::write.table(GenoOUT, file = OutFile,
                         row.names = TRUE, col.names = FALSE, quote = FALSE)
    }
  } else {  # output to console. invisible() to only print if assigned
    if (OutFormat == 'ped') {
      invisible( list(ped = cbind(FamOUT, GenoOUT), map = MapOUT) )
    } else {
      invisible( GenoOUT )
    }
  }
}



#######################################################################
#######################################################################

#' @title Extract Sex and Birth Year from PLINK File
#'
#' @description Convert the first six columns of a PLINK .fam, .ped or
#'  .raw file into a three-column lifehistory file for sequoia. Optionally
#'   FID and IID are combined.
#'
#' @details The first 6 columns of PLINK .fam, .ped and .raw files are by
#' default FID - IID - father ID (ignored) - mother ID (ignored) - sex -
#' phenotype.
#'
#' @param PlinkFile character string with name of genotype file to be converted.
#' @param UseFID use the family ID column. The resulting ids (rownames of GenoM)
#'   will be in the form FID__IID.
#' @param SwapSex change the coding from PLINK default (1=male, 2=female) to
#'   sequoia default (1=female, 2=male); any other numbers are set to NA.
#' @param FIDsep characters inbetween FID and IID in composite-ID. By default a
#'   double underscore is used, to avoid problems when some IIDs contain an
#'   underscore. Only used when UseFID=TRUE.
#' @param LifeHistData  dataframe with additional sex and birth year info. In
#'   case of conflicts, LifeHistData takes priority, with a warning. If
#'   UseFID=TRUE, IDs in LifeHistData are assumed to be already as FID__IID.
#'
#' @return A dataframe with id, sex and birth year, which can be used as input
#'  for \code{\link{sequoia}}.
#'
#' @seealso \code{\link{GenoConvert}}, \code{\link{PedStripFID}} to reverse
#'  \code{UseFID}.
#'
#' @examples
#' \dontrun{
#' # combine FID and IID in dataframe with additional sex & birth years
#' ExtraLH$FID_IID <- paste(ExtraLH$FID, ExtraLH$IID, sep = "__")
#' LH.new <- LHConvert(PlinkFile, UseFID = TRUE, FIDsep = "__",
#'                     LifeHistData = ExtraLH)
#' }
#'
#' @export

LHConvert <- function(PlinkFile = NULL, UseFID = FALSE,
                      SwapSex = TRUE, FIDsep="__", LifeHistData=NULL)
{
  if (is.null(PlinkFile)) stop("please provide 'InFile'")
  if (!file.exists(PlinkFile)) stop("cannot find 'PlinkFile'")
  if (UseFID & FIDsep %in% c("", " ", "\t", "\n")) stop("sep can not be whitespace")
  if (!is.null(LifeHistData)) {
    LHIN <- CheckLH(LifeHistData, sorted=FALSE)
  }

  ncol <- length(scan(PlinkFile, nlines=1, what="real", quiet=TRUE))
  TMP <- scan(PlinkFile, skip=1, what=as.list(c(rep("character", 2), rep("numeric", 4),
                                                rep("NULL", ncol-6))), quiet=TRUE)

  LH <- data.frame(id = TMP[[2]],
                   Sex = TMP[[5]],
                   BirthYear = TMP[[6]],
                   stringsAsFactors=FALSE)
  if (SwapSex) {
    LH$Sex <- ifelse(LH$Sex==1, 2,
                     ifelse(LH$Sex==2, 1,
                            NA))
  }

  if (UseFID) {
    IDX <- data.frame(id.old = TMP[[2]],
                      id.new = paste(TMP[[1]], TMP[[2]], sep=FIDsep),
                      stringsAsFactors=FALSE)
    LH <- merge(LH, IDX, by.x="id", by.y="id.old", all.x=TRUE)
    LH$id <- ifelse(!is.na(LH$id.new), LH$id.new, LH$id)
    LH <- LH[, c("id", "Sex", "BirthYear")]
  }

  if (!is.null(LHIN)) {
    names(LHIN) <- c("id", "Sex", "BirthYear")
    LH$Sex[!LH$Sex %in% c(1,2)] <- NA
    LHIN$Sex[!LHIN$Sex %in% c(1,2)] <- NA
    LH$BirthYear[LH$BirthYear < 0] <- NA
    LHIN$BirthYear[LHIN$BirthYear < 0] <- NA

    chk <- merge(LH, LHIN, by="id")
    n.sexmismatch <- sum(chk$Sex.x != chk$Sex.y, na.rm=T)
    n.BYmismatch <- sum(chk$BirthYear.x != chk$BirthYear.y, na.rm=T)
    if (n.sexmismatch > 0) {
      these <- with(chk, id[which(!is.na(Sex.x) & !is.na(Sex.y) & Sex.x!=Sex.y)])
      cli::cli_alert_warning("There are {n.sexmismatch} sex mismatches: ")
      cli::cli_li(these[1:min(length(these), 3)])
      if (length(these) > 3)  cli::cli_li("...")
    }
    if (n.BYmismatch > 0) {
      these <- with(chk, id[which(!is.na(BirthYear.x) & !is.na(BirthYear.y) & BirthYear.x!=BirthYear.y)])
      cli::cli_alert_warning("There are {n.BYmismatch} birth year mismatches: ")
      cli::cli_li(these[1:min(length(these), 3)])
      if (length(these) > 3)  cli::cli_li("...")
    }

    LH <- MergeFill(LH, LHIN, by="id", overwrite=TRUE, all=TRUE)
  }

  LH
}

