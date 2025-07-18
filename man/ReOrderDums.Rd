% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CalcParentProbs.R
\name{ReOrderDums}
\alias{ReOrderDums}
\title{Re-order likelihood vector}
\usage{
ReOrderDums(LLM, Complex = "full")
}
\arguments{
\item{LLM}{a matrix with log10-likelihoods; individuals in rows.}
}
\value{
a matrix with similar dimensions, but changed column order:
  \item{.}{GP \eqn{\rightarrow} PO (GP of sibship = parent of dummy)}
  \item{.}{FA \eqn{\rightarrow} FS}
  \item{.}{HA \eqn{\rightarrow} HS}
  \item{.}{FS \eqn{\rightarrow} FA (FS of sibs=offspring of dummy, NOT parent)}
  \item{.}{HS \eqn{\rightarrow} HA}
  \item{.}{U \eqn{\rightarrow} U}
  \item{.}{PO \eqn{\rightarrow} DUP (only if 'DUP' already among column names)}
}
\description{
change the order for dummy individuals so that the
 results are easily comparable with genotyped individuals
}
\examples{
# TODO

}
\keyword{internal}
