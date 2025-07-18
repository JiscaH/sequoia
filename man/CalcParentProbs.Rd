% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CalcParentProbs.R
\name{CalcParentProbs}
\alias{CalcParentProbs}
\title{Calculate assignment probabilities}
\usage{
CalcParentProbs(Pedigree = NULL, GenoM = NULL, quiet = FALSE, nCores = 1, ...)
}
\arguments{
\item{Pedigree}{dataframe with columns id-dam-sire. By default, any
non-genotyped individuals are 'dummified'; use \code{Module='par'} to
ignore them.}

\item{GenoM}{numeric matrix with genotype data: One row per individual,
one column per SNP, coded as 0, 1, 2, missing values as a negative number
or NA. You can reformat data with \code{\link{GenoConvert}}, or use other
packages to get it into a genlight object and then use \code{as.matrix}.}

\item{quiet}{logical, suppress messages. No progress is printed when >1 core
is used.}

\item{nCores}{number of computer cores to use. If \code{2} or \code{4},
package \pkg{parallel} is used (other values are not applicable).}

\item{...}{Additional arguments passed to \code{\link{CalcPairLL}}, such as
the genotyping error rate \code{Err}, age information in
\code{LifeHistData} and \code{AgePrior}, or \code{InclDup} to include the
probability that the two samples are duplicates.}
}
\value{
the \code{Pedigree} dataframe with the three applicable columns
renamed to id-dam-sire, and 7 additional columns:
\item{Probdam}{Probability that individual in dam column is the maternal parent,
  rather than otherwise related (LL(PO)/sum(LL))}
\item{Probsire}{Analogous for sire}
\item{Probpair}{Probability for id-dam-sire trio. Approximated as the minimum
  of dam conditional on sire and sire conditional on dam, thus not including
  e.g. both being siblings (those other configurations are considered by
  sequoia during pedigree reconstruction, but can (currently) not be accessed
  directly)}
\item{dam_alt, sire_alt}{Most likely alternative (not PO) relationship
 between id-dam and id-sire, respectively}
\item{Probdam_alt, Probsire_alt}{Probability of most likely alternative
 relationship}
}
\description{
For each assigned offspring-parent pair, calculate the
 probability they are parent-offspring vs otherwise related. Probabilities
 are scaled to sum to one across all possible* relationships between the pair
 or trio; see Details.
}
\details{
The returned probabilities are calculated from the likelihoods used
  throughout the rest of this package, by scaling them to sum to one across
  all possible relationships. For \code{Complex='simp'} these are
  PO=parent-offspring, FS=full siblings, HS=half siblings, GP=grand-parental,
  FA=full avuncular, HA=third degree relatives (incl half avuncular), and
  U=unrelated. For \code{Complex='full'} there are numerous double relationship
  considered (PO & HS, HS & HA, etc), making both numerator and denominator
  in the scaling step less unambiguous, and the returned probabilities an
  approximation.

 The likelihoods are calculated by calling \code{\link{CalcPairLL}} once or
 twice for each id-dam and id-sire pair: once not conditioning on the
 co-parent, and once conditional on the co-parent, if any. For genotyped
 individuals this is done with \code{focal='PO'}, and for dummy individuals
 with \code{focal='GP'}.

 For relationships between a genotyped and a dummy individual, it may only be
 possible to determine that the genotyped individual is a second degree
 relative (GP, HS, or FA) to the dummy's offspring. This then results in a
 probability of at most 0.33, even when the two are indeed parent and
 offspring.

 See \code{\link{CalcPairLL}} and the vignettes for further details.

Note that for large pedigrees this function can be fairly slow, especially
when using \code{\link{CalcPairLL}}'s default \code{Module='ped'} and
\code{Complex='full'}.

Subsetting the genotype data may give different results, as the likelihoods
 and thus the probabilities depend on the allele frequencies in the sample.
}
\examples{

test_ped <- Ped_griffin[21:25,]
# add an incorrect sire to illustrate
test_ped$sire <- as.character(test_ped$sire)
test_ped$sire[5] <- 'i057_2003_M'
Ped_with_probs <- CalcParentProbs(test_ped, Geno_griffin)
print(Ped_with_probs, digits=2)
# Any non-genotyped non-'dummifiable' individuals are automatically skipped

# To get likelihoods for 'all' relationships, not just probabilities for
# PO & (next-)most-likely:
LL_sire_single <- CalcPairLL(
  Pairs = data.frame(id1=test_ped$id,
                     id2=test_ped$sire,
                     dropPar1='both', # drop both -> id2 as single parent
                     focal='PO'),
  Pedigree = Ped_griffin,   # pedigree to condition on
  GenoM = Geno_griffin, Plot=FALSE)
}
\seealso{
\code{\link{CalcPairLL}}, \code{\link{LLtoProb}}
}
