#=======================================================================
#' @title Check LifeHistData
#'
#' @description Check that the provided LifeHistData is in the correct format.
#'
#' @param LifeHistData the dataframe with ID - Sex - Birth year, and optionally
#'   BY.min - BY.max - YearLast.
#' @param gID  character vector with names of genotyped individuals, i.e.
#'   rownames(GenoM).
#' @param sorted  logical, return lifehistdata for genotyped individuals only,
#'   in strictly the same order. Will including padding with 'empty' rows if an
#'   individual in gID was not in the input-LH.
#' @param returnDups  logical, instead of just the (sorted) LifeHistData, return
#'   a list that also includes a dataframe with duplicate entries and/or a
#'   character vector with genotyped IDs not occuring in LifeHistData.
#'
#' @return A dataframe with LifeHistData formatted for use by the Fortran
#'    part of the program, or a list with duplicate and missing entries.
#'
#' @keywords internal
#' @noRd

CheckLH <- function(LifeHistData, gID = NA, sorted=TRUE, returnDups = FALSE)
{
  if (is.null(gID))  gID <- NA
  if (is.null(LifeHistData)) {
    LH.OUT <- data.frame(ID = gID,
                         Sex = 3,
                         BirthYear = -999,
                         BY.min = -999,
                         BY.max = -999,
                         Year.last = -999,
                         stringsAsFactors = FALSE)
    if (returnDups) {
      return( list(DupLifeHistID = NULL, NoLH = NULL, LifeHistData = LH.OUT) )
    } else {
      return( LH.OUT )
    }
  }
  if (all(is.na(LifeHistData)))
    stop("invalid value for LifeHistData, provide NULL or dataframe",
         call.=FALSE)

  if (ncol(LifeHistData) < 3)
    stop("LifeHistData must have at least 3 columns: ID - Sex - BirthYear (or be NULL)")


  # check IDs (column 1) ---
  colnames(LifeHistData)[1] <- 'ID'
  LifeHistData$ID <- as.character(LifeHistData$ID)
  if (any(grepl(" ", LifeHistData$ID))) {
    NotOK <- LifeHistData$ID[grepl(" ", LifeHistData$ID)]
    cli::cli_alert_danger(c("LifeHistData IDs (column 1) contains {length(NotOK)}",
                                " IDs that include spaces:"))
    cli::cli_li(paste0("'", NotOK[1:min(length(NotOK), 3)], "'"))   # quotation marks to see leading/trailing blanks
    if (length(NotOK) > 3)  cli::cli_li("...")
    stop("IDs may only contain alphanumerical characters and underscore", call.=FALSE)
  }
  if (!all(is.na(gID))) {
    gID <- as.character(gID)
    if (length(intersect(LifeHistData$ID, gID))==0) {
      cli::cli_alert_danger("Incorrect LifeHistData object, or incorrect format")
      stop("None of the IDs in LifeHistData column 1 match the rownames of GenoM",
           call.=FALSE)
    }
  } else {
    if (sorted)  stop("if sorted=TRUE, gID (ordered IDs in GenoM) must be provided", call.=FALSE)
  }


  # Flexible column order and capitalisation ---

  which_else <- function(x, z) ifelse(any(x), which(x), z)   # z = default column if name not found

  LH_colnames <- tolower(colnames(LifeHistData))
  colnum <- c('Sex' = which_else(grepl('sex', LH_colnames), 2),
              'BirthYear' = which_else((grepl('birthyear', LH_colnames) | grepl('by', LH_colnames)) &
                                         !grepl('min', LH_colnames) & !grepl('max', LH_colnames), 3),
              'BY.min' = which_else((grepl('birthyear', LH_colnames) | grepl('by', LH_colnames)) &
                                      grepl('min', LH_colnames), 4),
              'BY.max' = which_else((grepl('birthyear', LH_colnames) | grepl('by', LH_colnames)) &
                                      grepl('max', LH_colnames), 5),
              'Year.last' = which_else(grepl('last', LH_colnames), 6) )
  if (any(duplicated(colnum))) {
    cli::cli_alert_danger("Column order and/or column names in LifeHistData is unclear")
    stop("Please provide LifeHistData as ",
         "ID - Sex - BirthYear (- BY.min - BY.max - Year.last)", call.=FALSE)
  }

  # optional columns
  for (x in c('BY.min', 'BY.max', 'Year.last')) {
    if (colnum[x] > ncol(LifeHistData)) {
      LifeHistData[, x] <- -999
    }
  }

  # column renaming & order ---
  for (x in seq_along(colnum)) {
    colnames(LifeHistData)[ colnum[x] ] <- names(colnum)[x]
  }
  LifeHistData <- LifeHistData[, c('ID', 'Sex', 'BirthYear', 'BY.min', 'BY.max', 'Year.last')]


  # check Sex ----
  if (!all(LifeHistData$Sex %in% 1:4)) {
    prop_sex_OK <- sum(LifeHistData$Sex %in% c(1:4,NA)) / nrow(LifeHistData)
    sex_msg <- "LifeHistData column 'Sex' should be coded as 1=female, 2=male, 3/<NA>=unknown, 4=hermaphrodite"
    if (prop_sex_OK < 0.1) {
      # recode from F/M/.. to 1/2/3?
      if (sum(LifeHistData$Sex %in% c('f','F','m','M'))/nrow(LifeHistData) > 0.1) {
        cli::cli_alert_warning("Recoding `LifeHistData` column `Sex` from F/M/.. to 1=female, 2=male, 3=unknown")
        LifeHistData$Sex[LifeHistData$Sex %in% c('f','F')] <- 1
        LifeHistData$Sex[LifeHistData$Sex %in% c('m','M')] <- 2
      } else {
        cli::cli_alert_danger(c('{(round((1-prop_sex_OK)*100)}% of entries',
        " in `LifeHistData` column `Sex` have an unrecognised coding"))
        stop(sex_msg, call.=FALSE)
      }
    } else {
      if (!all(LifeHistData$Sex %in% c(1:4,NA))) {
        cli::cli_alert_warning(c(sex_msg, "\nThe following values are converted to 3=unknown:"))
        cli::cli_li(setdiff(LifeHistData$Sex, c(1:4,NA)))
      }
    }
    LifeHistData$Sex[!LifeHistData$Sex %in% c(1:4)] <- 3
  }


  # check birth years ----
  LifeHistData$BY.max[which(LifeHistData$BY.max < 0)] <- NA
  if (any(LifeHistData$BY.max < LifeHistData$BY.min, na.rm=TRUE)) {
    stop("'BY.max' must be greater than or equal to 'BY.min', or <NA>/-999", call.=FALSE)
  }

  for (x in c("BirthYear", "BY.min", "BY.max", 'Year.last')) {
    IsInt <- check.integer(LifeHistData[,x])
    if (any(!IsInt, na.rm=TRUE)) {
      cli::cli_alert_warning("In `LifeHistData` column `{x}`, the following values are converted to <NA>/-999: ")
      cli::cli_li(unique(LifeHistData[,x][IsInt %in% FALSE]))
    }
    LifeHistData[, x] <- ifelse(IsInt, suppressWarnings(as.integer(as.character(LifeHistData[, x]))), NA)
    LifeHistData[is.na(LifeHistData[,x]), x] <- -999
    LifeHistData[LifeHistData[,x] < 0, x] <- -999
  }


  # duplicate check ----
  LifeHistData <- unique(LifeHistData)    # doesn't matter if both entries are identical
  OUT <- list(DupLifeHistID = NULL, NoLH = NULL)
  if (any(duplicated(LifeHistData$ID))) {
    if (returnDups) {
      r1 <- which(duplicated(LifeHistData$ID, fromLast=TRUE))
      r2 <- which(duplicated(LifeHistData$ID))
      LHdup <- merge(LifeHistData[r1,], LifeHistData[r2,], by='ID', suffixes=c('.1', '.2'))
      col_order <- c(t(outer(c('Sex', 'BirthYear', 'BY.min', 'BY.max', 'Year.last'), 1:2, paste, sep='.')))
      OUT$DupLifeHistID <- data.frame(ID = LHdup$ID,
                                      row1 = r1,
                                      row2 = r2,
                                      LHdup[, c('ID', intersect(col_order, colnames(LHdup)))],
                                      stringsAsFactors = FALSE)
    }
    cli::cli_alert_warning("duplicate IDs found in `LifeHistData`, first entry will be used")
    LifeHistData <- LifeHistData[!duplicated(LifeHistData$ID), ]
  }
  if (returnDups & !all(is.na(gID))) {
    NoLH <- setdiff(gID, LifeHistData$ID)
    if (length(NoLH)>0)  OUT$NoLH <- NoLH
  }

  # return results ----
  if (sorted)
    LifeHistData <- orderLH(LifeHistData, gID)

  if (returnDups) {
    return( c(OUT, list(LifeHistData=LifeHistData)) )
  } else {
    return( LifeHistData )
  }
}


#=======================================================================
#=======================================================================
#' @title Order Lifehistory Data
#'
#' @description Order lifehistory data to match order of IDs in genotype data,
#'   filling in gaps with missing values.
#'
#' @param LH dataframe with lifehistory information:
#' \describe{
#'  \item{ID}{max. 30 characters long,}
#'  \item{Sex}{1 = females, 2 = males, other numbers = unknown,}
#'  \item{Birth Year}{(or hatching year) Use negative numbers to denote
#'  missing values.}
#'  \item{BY.min}{minimum birth year (optional)}
#'  \item{BY.max}{maximum birth year (optional)}}
#' @param gID character vector with IDs in genotype data, in order of
#'   occurrence.
#'
#' @return A dataframe with the same 5 columns, but with individuals in exactly
#'   the same order as gID, including padding with 'empty' rows if an individual
#'   in gID was not in the input-LH. Missing values are recoded to 3 for the
#'   'Sex' column, and -999 for the birth year columns.
#'
#' @keywords internal 
#' @noRd

orderLH <- function(LH=NULL, gID=NULL) {
  if (!all(gID %in% LH$ID)) {
    LHX <- data.frame(ID = gID[!gID %in% LH$ID],
                      Sex = 3,
                      BirthYear = -999,
                      BY.min = -999,
                      BY.max = -999,
                      Year.last = -999)
  }
  if (is.null(LH)) {
    LHF <- LHX
  } else {
    if (!all(gID %in% LH$ID)) {
      LH <- merge(LHX, LH, all = TRUE)
    }
    LHF <- LH[match(gID, LH$ID), ]
    # fill in gaps: some gID may not be in LH
    LHF$Sex[is.na(LHF$Sex)] <- 3
    for (y in c("BirthYear", "BY.min", "BY.max", "Year.last")) {
      LHF[is.na(LHF[,y]), y] <- -999
    }
  }
  return( LHF )
}
