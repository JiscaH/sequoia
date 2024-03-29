% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SeqListSummary.R
\name{PlotSeqSum}
\alias{PlotSeqSum}
\title{Plot Summary Overview of sequoia Output}
\usage{
PlotSeqSum(SeqSum, Pedigree = NULL, Panels = "all", ask = TRUE)
}
\arguments{
\item{SeqSum}{list output from \code{\link{SummarySeq}}.}

\item{Pedigree}{dataframe with at least id, dam and sire in columns 1-3,
respectively. If columns with parental LLRs and/or Mendelian errors are
present, these will be plotted as well.}

\item{Panels}{character vector with panel(s) to plot. Choose from 'all',
'G.parents' (parents of genotyped individuals), 'D.parents' (parents of
dummies), 'O.parents' (parents of non-genotyped non-dummies), sibships',
'LLR', 'OH'.}

\item{ask}{ask for user key stroke before proceeding to next plot.}
}
\description{
visualise the numbers of assigned parents, sibship sizes, and
  parental LLRs
}
\examples{
sumry <- SummarySeq(SeqOUT_griffin, Plot=FALSE)
PlotSeqSum(sumry, SeqOUT_griffin$Pedigree, Panels='all', ask=FALSE)

}
