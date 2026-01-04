#' @title Write Sequoia Output to File
#'
#' @description The various list elements returned by \code{sequoia} are each
#'   written to text files in the specified folder, or to separate sheets in a
#'   single excel file (requires library \pkg{openxlsx}).
#'
#' @details The text files can be used as input for the stand-alone Fortran
#'   version of sequoia, e.g. when the genotype data is too large for R. See
#'   \code{vignette('sequoia')} for further details.
#'
#' @param SeqList list returned by \code{\link{sequoia}}, to be written out.
#' @param GenoM  matrix with genetic data (optional). Ignored if
#'   OutFormat='xls', as the resulting file could become too large for excel.
#' @param MaybeRel list with results from \code{\link{GetMaybeRel}} (optional).
#' @param PedComp list with results from \code{\link{PedCompare}} (optional).
#'   \code{SeqList$DummyIDs} is combined with \code{PedComp$DummyMatch} if both
#'   are provided.
#' @param OutFormat 'xls' or 'txt'.
#' @param folder the directory where the text files will be written; will be
#'   created if it does not already exists. Relative to the current working
#'   directory, or NULL for current working directory. Ignored if
#'   \code{OutFormat='xls'}.
#' @param file the name of the excel file to write to, ignored if
#'   \code{OutFormat='txt'}.
#' @param quiet suppress messages.
#'
#' @seealso \code{\link{writeColumns}} to write to a text file, using white
#'   space padding to keep columns aligned.
#'
#' @examples
#' \dontrun{
#' writeSeq(SeqList, OutFormat="xls", file="MyFile.xlsx")
#'
#' # add additional sheet to the excel file:
#' library(openxlsx)
#' wb <- loadWorkbook("MyFile.xlsx")
#' addWorksheet(wb, sheetName = "ExtraData")
#' writeData(wb, sheet = "ExtraData", MyData, rowNames=FALSE)
#' saveWorkbook(wb, "MyFile.xlsx", overwrite=TRUE, returnValue=TRUE)
#'
#' # or: (package requires java & is trickier to install)
#' xlsx::write.xlsx(MyData, file = "MyFile.xlsx", sheetName="ExtraData",
#'       col.names=TRUE, row.names=FALSE, append=TRUE, showNA=FALSE)
#' }
#'
#' @export

writeSeq <- function(SeqList,
                     GenoM = NULL,
                     MaybeRel = NULL,
                     PedComp = NULL,
                     OutFormat = "txt",
                     folder = "Sequoia-OUT",
                     file = "Sequoia-OUT.xlsx",
                     quiet = FALSE) 
{
  if (!(isTRUE(quiet) | isFALSE(quiet)))  stop("'quiet' must be TRUE or FALSE")
  if (!OutFormat %in% c("xls", "xlsx", "txt"))  stop("Invalid OutFormat")
  if (!inherits(SeqList, 'list'))  stop("SeqList should be a list")

  if (!is.null(MaybeRel)) {
    for (maybe in c("MaybeRel", "MaybeParent", "MaybeTrio")) {
      if (maybe %in% names(MaybeRel)) {
        SeqList[[maybe]] <- MaybeRel[[maybe]]
      }
    }
  }

  if (!is.null(GenoM)) {
    if (!is.matrix(GenoM) & !is.data.frame(GenoM))  stop("GenoM should be a matrix")
    if ("PedigreePar" %in% names(SeqList)) {
      if (!all(rownames(GenoM) %in% SeqList$PedigreePar$id)) {
        stop("Not all ids in 'GenoM' occur in SeqList$PedigreePar")
      } else if (!all(rownames(GenoM) == SeqList$PedigreePar$id)) {
        stop("rownames of 'GenoM' do not match order of ids in SeqList$PedigreePar")
      }
    }
    if ("LifeHist" %in% names(SeqList)) {
      if (length(intersect(SeqList[["LifeHist"]]$ID, rownames(GenoM))) == 0) {
        stop("rownames of 'GenoM' shares no common IDs with 'LifeHist' in SeqList")
      }
    }
    if(nrow(GenoM)!= SeqList$Specs$NumberIndivGenotyped) {
      ANS <- readline(prompt = paste("Number of individuals according to Specs differs",
      " from number of rows in GenoM (", SeqList$Specs$NumberIndivGenotyped, "/", nrow(GenoM),
      ").\n Press Y to continue and fix manually in `SequoiaSpecs.txt' "))
      if (!substr(ANS, 1, 1) %in% c("Y", "y")) {
        stop(call.=FALSE)
      }
    }
    if(ncol(GenoM)!= SeqList$Specs$NumberSnps) {
      stop(paste("Number of SNPs according to Specs differs from number of rows in GenoM (",
                 SeqList$Specs$NumberSnps, "vs", ncol(GenoM), ")"))
    }
  } else if (OutFormat == "txt") {
    cli::cli_alert_warning("No `GenoM` specified")
  }

  # write excel file ----
  if (OutFormat == "xlsx") OutFormat <- "xls"
  if (grepl('xls', file) & !grepl('Sequoia-OUT', file))  OutFormat <- "xls"
  if (OutFormat == "xls") {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      if (interactive() & !quiet) {
        ANS <- readline(prompt = paste("library 'openxlsx' not found. Install Y/N? "))
        if (!substr(ANS, 1, 1) %in% c("Y", "y")) stop(call.=FALSE)
      }
      utils::install.packages("openxlsx")
    }
    write.seq.xls(SeqList, file=file, PedComp=PedComp, quiet=quiet)


  # txt: check if folder & files exist ----
  } else if (OutFormat == "txt") {
  curdir <- getwd()
  if (is.null(folder)) folder = curdir
  dir.create(folder, showWarnings = FALSE)
  setwd(folder)
  if (any(file.exists("Geno.txt", "Specs.txt", "AgePriors.txt", "Parents.txt",
                      "Readme.txt")) &
      interactive() & !quiet) {
    ANS <- readline(prompt = paste("Writing data to '", folder,
                                   "' will overwrite existing file(s). Continue Y/N? "))
    if (!substr(ANS, 1, 1) %in% c("Y", "y", "J", "j", "")) {
      setwd(curdir)
      stop(call.=FALSE)
    }
  }

  # prep specs ----
  SpecsOUT <- cbind("Genofile" = "Geno.txt",
                    "LHfile" = "LifeHist.txt",
                    "nIndLH" = nrow(SeqList[["LifeHist"]]),
                    SeqList$Specs)
  SpecsOUT$Complexity <- c("mono"=0, "simp"=1, "full"=2)[SpecsOUT$Complexity]
  SpecsOUT$Herm <- c("no"=0, "A"=1, "B"=2)[SpecsOUT$Herm]
  SpecsOUT$UseAge <- c("extra"=2, "yes"=1, "no"=0)[SpecsOUT$UseAge]
  if ('FindMaybeRel' %in% names(SpecsOUT)) {   # dropped from version 2.4
    SpecsOUT$FindMaybeRel <- as.numeric(SpecsOUT$FindMaybeRel)
  }
  SpecsOUT$CalcLLR <- as.numeric(SpecsOUT$CalcLLR)
  for (x in c("SequoiaVersion", "TimeStart", "TimeEnd")) {
    if (!x %in% names(SpecsOUT))  next   # if SeqList from version 1.x
    SpecsOUT[[x]] <- as.character(SpecsOUT[[x]])
  }
  SpecsOUT <- SpecsOUT[!names(SpecsOUT) %in% "Module"]  # run-time option for stand-alone sequoia.


  # write text files ----
  OPT <- options()
  options(scipen = 10)
  write(paste("These files were created by R package sequoia on ", date()),
        file="README.txt")
  if (!is.null(GenoM)) utils::write.table(GenoM, file="Geno.txt",
              quote=FALSE, row.names=TRUE, col.names=FALSE)
  utils::write.table(as.data.frame(t(SpecsOUT)), file="SequoiaSpecs.txt",
              sep = "\t,\t", quote = FALSE, col.names = FALSE)
  writeColumns(SeqList$AgePriors, "AgePriors.txt", row.names=FALSE)
  if (ncol(SeqList[["LifeHist"]])>6) {
    writeColumns(SeqList[["LifeHist"]][,1:6], "LifeHist.txt", row.names=FALSE)
  } else {
    writeColumns(SeqList[["LifeHist"]], "LifeHist.txt", row.names=FALSE)
  }

  if ("PedigreePar" %in% names(SeqList)) {
    write.parents(ParentDF = SeqList$PedigreePar, LifeHistData = SeqList[["LifeHist"]],
                  GenoM = GenoM)
		# compatable to be read in by stand-alone sequoia
  }

  if ("Pedigree" %in% names(SeqList)) {
    writeColumns(SeqList[["Pedigree"]], "Pedigree.txt", row.names=FALSE)
  }

  if ("DummyIDs" %in% names(SeqList)) {
    if ("DummyMatch" %in% names(PedComp)) {
      Dummies <- merge(PedComp$DummyMatch, SeqList$DummyIDs)
    } else {
      Dummies <- SeqList$DummyIDs
    }
    writeColumns(Dummies, "DummyIDs.txt", row.names=FALSE)
  }

  if ("MaybeRel" %in% names(SeqList)) {
    MaybePairs <- SeqList$MaybeRel
  } else if ("MaybeParent" %in% names(SeqList)) {
    MaybePairs <- SeqList$MaybeParent
  } else {
    MaybePairs <- NULL
  }
  if (!is.null(MaybePairs)) {
    writeColumns(MaybePairs, "MaybePairs.txt", row.names=FALSE)
  }
  if (!is.null(SeqList$MaybeTrio)) {
    writeColumns(SeqList$MaybeTrio, "MaybeTrios.txt", row.names=FALSE)
  }

  if (!is.null(PedComp)) {
    Counts.out <- cbind(Parent = rep(c("dam", "sire"), each=7),
                        Cat = rep(dimnames(PedComp$Counts)[[1]], 2))
    Counts.out <- cbind(Counts.out, rbind(PedComp$Counts[,,"dam"],
                                      PedComp$Counts[,,"sire"]))
    writeColumns(Counts.out, "PedCompare-Counts.txt", row.names=FALSE)

    writeColumns(PedComp$MergedPed, "PedCompare-MergedPed.txt",
                     row.names=FALSE)
  }

  options(OPT)
  setwd(curdir)
  if(!quiet) cli::cli_alert_info("Output written to {normalizePath(folder, winslash='/')}")
  } else {
    stop("OutFormat not supported.")
  }
}


#==========================================================================
write.parents <- function(ParentDF, LifeHistData, GenoM=NULL, file="Parents.txt") {
  names(ParentDF)[1:3] <- c("id", "dam", "sire")
  names(LifeHistData)[1:3] <- c("ID", "Sex", "BirthYear")

  if (!is.null(GenoM)) {
    gID <- rownames(GenoM)
  } else {
    gID <- ParentDF$id
  }

  Par <- MergeFill(ParentDF,
                   data.frame(id = gID,
                              LLRdam = NA, LLRsire = NA, LLRpair = NA,
                              OHdam = NA, OHsire = NA, MEpair = NA,
                              rowid = NA, rowdam = NA, rowsire = NA,
                              stringsAsFactors=FALSE),
                   by = "id", all = TRUE)
  rownames(Par) <- as.character(Par$id)
  Par <- Par[gID, ]    # merge ignores sort=FALSE
  for (x in c("dam", "sire")) Par[is.na(Par[, x]), x] <- "NA"
  for (x in c("LLRdam", "LLRsire", "LLRpair")) {
    Par[is.na(Par[, x]), x] <- 999
  }
  for (x in c("OHdam", "OHsire", "MEpair")) {
    Par[is.na(Par[, x]), x] <- -9
  }
  for (x in c("id", "dam", "sire")) {
    Par[, paste0("row", x)] <- as.numeric(factor(Par[,x], levels=gID))
    Par[is.na(Par[, paste0("row", x)]), paste0("row", x)] <- 0
  }
  if (any(is.na(Par$id))) {
    if (!is.null(GenoM)) {
      stop("Some id's in PedigreePar do not occur in GenoM!")
    } else {
      stop("Some id's in PedigreePar are <NA>")
    }
  }
  writeColumns(Par, file = file, row.names=FALSE)
}



#==========================================================================


write.seq.xls <- function(SeqList, file, PedComp=NULL, quiet) {

  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("library 'openxlsx' not found")
  }

  # check if file exists
  if (file.exists(file) & interactive() & !quiet) {
    ANS <- readline(prompt = paste("Writing data to '", file,
                                   "' will overwrite existing file. Continue Y/N? "))
    if (!substr(ANS, 1, 1) %in% c("Y", "y")) stop(call.=FALSE)
  }

  out_list <- list()   # do in wb straight away (?)
  if ("Pedigree" %in% names(SeqList)) {
    out_list[['Pedigree']] <- SeqList[["Pedigree"]]
  }
  if ('PedigreePar' %in% names(SeqList)) {
    out_list[['Parents']] <- SeqList[['PedigreePar']]
  }
  if (!any(c("Pedigree", "PedigreePar") %in% names(SeqList))) {
    out_list[['Parents']] <- 'None'
  }
  out_list[['Run parameters']] <- t(SeqList$Specs)

  Names <- c("LifeHist", "DupGenotype", "DupLifeHistID",
             "AgePriors", "DummyIDs", "TotLikPar", "TotLikSib",
             "MaybeRel", "MaybeParent", "MaybeTrio")
  for (x in Names) {
    if (!x %in% names(SeqList))  next
    if (is.null(dim(SeqList[[x]])))  next
    if (nrow(SeqList[[x]])==0)  next
    if (x == "DummyIDs") {
      if ("DummyMatch" %in% names(PedComp)) {
        SeqList$DummyIDs <- merge(PedComp$DummyMatch, SeqList$DummyIDs)
      }
    }
    out_list[[x]] <- SeqList[[x]]

  }

  if (!is.null(PedComp)) {  # pedcompare output
    Counts.out <- cbind(Parent = rep(c("dam", "sire"), each=7),
                        Cat = rep(dimnames(PedComp$Counts)[[1]], 2))
    Counts.out <- cbind(Counts.out, rbind(PedComp$Counts[,,"dam"],
                                          PedComp$Counts[,,"sire"]))
    rownames(Counts.out) <- 1:nrow(Counts.out)  # not used, but otherwise warning

    out_list[['PedCompare-Counts']] <- Counts.out
    out_list[['PedCompare-MergedPed']] <- PedComp$MergedPed
  }

  ShowRowName <- rep(FALSE, length(out_list))
  names(ShowRowName) <- names(out_list)
  ShowRowName[c('Run parameters', 'AgePriors')] <- TRUE

  # openxlsx::write.xlsx does not give error if xlsx exists & is open
  # --> turn into workbook & use saveWorkbook instead
  wb <- openxlsx::buildWorkbook(out_list,
                                colNames = TRUE, na.string = '',
                                rowNames = ShowRowName,
                                headerStyle = openxlsx::createStyle(textDecoration = 'bold'))
  WriteSuccess <- openxlsx::saveWorkbook(wb, file,
                               overwrite=TRUE, # overwrite warning generated by writeSeq()
                               returnValue=TRUE)

  if (isTRUE(WriteSuccess)) {
    if(!quiet) cli::cli_alert_info("Output written to {file}")
  } else {
    stop(WriteSuccess$message, call.=FALSE)
  }
}



#======================================================================
#' @title Write Data to a File Column-wise
#'
#' @description Write data.frame or matrix to a text file, using white
#' space padding to keep columns aligned as in \code{print}.
#'
#' @param x the object to be written, preferably a matrix or data frame.
#'  If not, it is attempted to coerce x to a matrix.
#' @param file a character string naming a file.
#' @param row.names a logical value indicating whether the row names of x are
#'  to be written along with x.
#' @param col.names a logical value indicating whether the column names of x
#'  are to be written along with x.
#'
#' @export

writeColumns <- function(x, file="", row.names=TRUE,
                         col.names=TRUE) {
  M <- as.matrix(x)
  if(col.names)  M <- rbind(colnames(x), M)
  if (row.names) {
    if (col.names) M <- cbind(c("", rownames(x)), M)
    if (!col.names) M <- cbind(rownames(x), M)
  }
  write(t(format(M)), ncolumns = ncol(M), file = file)
}
