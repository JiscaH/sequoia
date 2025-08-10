#' @title Count Generations
#'
#' @description For each individual in a pedigree, count the number of
#'   generations since its most distant pedigree founder.
#'
#' @param Ped  dataframe, pedigree with the first three columns being id - dam -
#'   sire. Column names are ignored, as are additional columns.
#' @param StopIfInvalid if a pedigree loop is detected, stop with an error
#'   (TRUE, default) or return the Pedigree, to see where the problem(s) occur.
#'
#' @return  A vector with the generation number for each individual, starting at
#'   0 for founders. Offspring of G0 X G0 are G1, offspring of G0 X G1 or G1 x
#'   G1 are G2, etc. \code{NA} indicates a pedigree loop where an individual is
#'   its own ancestor (or that the pedigree has >1000 generations).
#'
#'   If no output name is specified, no results are returned, only an error
#'   message when the pedigree contains a loop.
#'
#'   To get more details about a pedigree loop, you can use
#'   https://github.com/JiscaH/sequoiaExtra/blob/main/find_pedigree_loop.R
#'
#' @seealso \code{\link{GetAncestors}, \link{GetDescendants}} to get all
#'   ancestors resp. descendants of a specific individual (with a warning if it
#'   is its own ancestor); \link{FindFamilies} to find connected sub-pedigrees.
#'
#' @examples
#' # returns nothing if OK, else error:
#' getGenerations(SeqOUT_griffin$Pedigree)
#'
#' # returns vector with generation numbers:
#' G <- getGenerations(SeqOUT_griffin$Pedigree, StopIfInvalid=FALSE)
#' table(G, useNA='ifany')
#' Ped_plus_G <- cbind(SeqOUT_griffin$Pedigree, G)
#'
#' @export

getGenerations <- function(Ped, StopIfInvalid=TRUE) {
  for (p in 2:3) {
    Ped[which(Ped[,p]==0), p] <- NA
  }
  Ped <- as.data.frame(Ped)

  if (nrow(Ped) == 0)  return( 0 )

  if (!all(na.exclude(unlist(Ped[,2:3])) %in% Ped[,1]))
    stop("Some parents do not occur in first column; please call PedPolish() first")

  Ped$gen <- NA  # individual's generation
  Ped$gen.dam <- NA
  Ped$gen.sire <- NA
  Ped$gen[is.na(Ped[,2]) & is.na(Ped[,3])] <- 0
  for (x in 0:1000) {
    Ped$gen.dam[is.na(Ped$gen.dam) & Ped[,2] %in% Ped[which(Ped$gen<=x), 1]] <- x
    Ped$gen.sire[is.na(Ped$gen.sire) & Ped[,3] %in% Ped[which(Ped$gen<=x), 1]] <- x
    Ped$gen[which(is.na(Ped$gen) &
                   (Ped$gen.dam<=x | is.na(Ped[,2])) &
                   (Ped$gen.sire<=x | is.na(Ped[,3])))] <- x+1
    if (!any(is.na(Ped$gen)))  break
  }
  if (any(is.na(Ped$gen))) {
    msg <- "An individual is its own ancestor, or >1000 generations."
    if (StopIfInvalid) {
      stop(msg,
           "\n Use getGenerations(, StopIfInvalid=FALSE) to find the culprit",
           call.=FALSE)
    } else {
      cli::cli_alert_danger(msg)
      cli::cli_alert_info("individuals within a pedigree loop have generation=<NA>")
    }
  }
  invisible( setNames(Ped$gen, Ped[,1]) )
}
