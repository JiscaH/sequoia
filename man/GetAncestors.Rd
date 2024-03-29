% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FindFamilies.R
\name{GetAncestors}
\alias{GetAncestors}
\title{Get ancestors}
\usage{
GetAncestors(id, Pedigree)
}
\arguments{
\item{id}{id of the individual}

\item{Pedigree}{dataframe with columns id - parent1 - parent2; only the
first 3 columns will be used.}
}
\value{
a list with as first element \code{id}, second parents, third
 grandparents, etc.. Each element is a vector with ids, the first three
 elements are named, the rest numbered. Ancestors are unsorted within each
 list element.
}
\description{
get all ancestors of an individual
}
\examples{
Anc_i200  <- GetAncestors('i200_2010_F', Ped_griffin)


}
