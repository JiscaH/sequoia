% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimGeno.R
\name{SimGeno}
\alias{SimGeno}
\title{Simulate Genotypes}
\usage{
SimGeno(
  Pedigree,
  nSnp = 400,
  ParMis = c(0, 0),
  MAF = 0.3,
  CallRate = 0.99,
  SnpError = 5e-04,
  ErrorFV = function(E) c((E/2)^2, E - (E/2)^2, E/2),
  ErrorFM = NULL,
  ReturnStats = FALSE,
  quiet = FALSE
)
}
\arguments{
\item{Pedigree}{dataframe, pedigree with the first three columns being id -
dam - sire, additional columns are ignored.}

\item{nSnp}{number of SNPs to simulate.}

\item{ParMis}{single number or vector length two with proportion of parents
with fully missing genotype. Ignored if CallRate is a named vector. NOTE:
default changed from 0.4 (up to version 2.8.5) to 0 (from version 2.9).}

\item{MAF}{either a single number with minimum minor allele frequency, and
allele frequencies will be sampled uniformly between this minimum
('min_maf') and 1-min_maf (NOTE: between min_maf and 0.5 up to version
2.11) OR a vector with minor allele frequency at each locus. In both cases,
this is the MAF among pedigree founders; the MAF in the sample will deviate
due to drift. Occasionally, alleles may get fixed, and all founders
carrying the allele simulated as non-genotyped, resulting in a sample MAF
of 0.}

\item{CallRate}{either a single number for the mean call rate (genotyping
success), OR a vector with the call rate at each SNP, OR a named vector
with the call rate for each individual. In the third case, ParMis is
ignored, and individuals in the pedigree (as id or as parent) not included
in this vector are presumed non-genotyped.}

\item{SnpError}{either a single value which will be combined with
\code{ErrorFV}, or a length 3 vector with probabilities (observed given
actual) hom|other hom, het|hom, and hom|het; OR a vector or 3XnSnp matrix
with the genotyping error rate(s) for each SNP.}

\item{ErrorFV}{function taking the error rate (scalar) as argument and
returning a length 3 vector with hom->other hom, hom->het, het->hom. May be
an 'ErrFlavour', e.g. 'version2.9'.}

\item{ErrorFM}{function taking the error rate (scalar) as argument and
returning a 3x3 matrix with probabilities that actual genotype i (rows) is
observed as genotype j (columns). See below for details. To use, set
\code{ErrorFV = NULL}}

\item{ReturnStats}{in addition to the genotype matrix, return the input
parameters and mean & quantiles of MAF, error rate and call rates.}

\item{quiet}{suppress messages.}
}
\value{
If \code{ReturnStats=FALSE} (the default), a matrix with genotype
  data in sequoia's input format, encoded as 0/1/2/-9.

  If \code{ReturnStats=TRUE}, a named list with three elements: list
  'ParamsIN', matrix 'SGeno', and list 'StatsOUT':
  \item{AF}{Frequency in 'observed' genotypes of '1' allele}
  \item{AF.act}{Allele frequency in 'actual' (without genotyping errors &
    missingness)}
  \item{SnpError}{Error rate per SNP (actual /= observed AND observed /=
    missing)}
  \item{SnpCallRate}{Non-missing per SNP}
  \item{IndivError}{Error rate per individual}
  \item{IndivCallRate}{Non-missing per individual}
}
\description{
Simulate SNP genotype data from a pedigree, with optional
  missingness, genotyping errors, and non-genotyped parents.
}
\details{
For founders, i.e. individuals with no
  known parents, genotypes are drawn according to the provided MAF and
  assuming Hardy-Weinberg equilibrium. Offspring genotypes are generated
  following Mendelian inheritance, assuming all loci are completely
  independent. Individuals with one known parent are allowed: at each locus,
  one allele is inherited from the known parent, and the other drawn from the
  genepool according to the provided MAF.
}
\section{Genotyping errors}{

  If \code{SnpError} is a length 3 vector, genotyping errors are generated
  following a length 3 vector with probabilities that 1) an actual homozygote
  is observed as the other homozygote, 2) an actual homozygote is observed as
  a heterozygote, and 3) an heterozygote is observed as an homozygote. The
  only assumption made is that the two alleles can be treated equally, i.e.
  observing actual allele $A$ as $a$ is as likely as observing actual $a$ as
  $A$.

  If \code{SnpError} is a single value, by default this is interpreted as a
  locus-level error rate (rather than allele-level), and equals the
  probability that a homozygote is observed as heterozygote, and the
  probability that a heterozygote is observed as either homozygote (i.e., the
  probability that it is observed as AA = probability that observed as aa =
  \code{SnpError}/2). The probability that one homozygote is observed as the
  other is (\code{SnpError}/2\eqn{)^2}. How this single value is rendered
  into a 3x3 error matrix is fully flexible and specified via \code{ErrorFM};
  see \code{link{ErrToM}} for details.

  The default values of \code{SnpError=5e-4} and \code{ErrorFM='version2.9'}
  correspond to the length 3 vector \code{c((5e-4/2)^2, 5e-4*(1-5e-4/2),
  5e-4/2)}.

  A beta-distribution is used to simulate variation in the error rate between
  SNPs, the shape parameter of this distribution can be specified via
  \code{\link{MkGenoErrors}}. It is also possible to specify the error rate
  per SNP.
}

\section{Call Rate}{

  Variation in call rates across SNPs is assumed to follow a highly skewed
  (beta) distribution, with many SNPs having call rates close to 1, and a
  narrowing tail of lower call rates. The first shape parameter defaults to 1
  (but see \code{\link{MkGenoErrors}}), and the second shape parameter is
  defined via the mean as \code{CallRate}. For 99.9\% of SNPs to have a call
  rate of 0.8 (0.9; 0.95) or higher, use a mean call rate of 0.969 (0.985;
  0.993).

  Variation in call rate between samples can be specified by providing a
  named vector to \code{CallRate}. Otherwise, variation in call rate and
  error rate between samples occurs only as side-effect of the random nature
  of which individuals are hit by per-SNP errors and drop-outs. Finer control
  is possible by first generating an error-free genotype matrix, and then
  calling \code{\link{MkGenoErrors}} directly on (subsets of) the matrix.
}

\section{Disclaimer}{
 This simulation is highly simplistic and assumes that
  all SNPs segregate completely independently, that the SNPs are in
  Hardy-Weinberg equilibrium in the pedigree founders. It assumes that
  genotyping errors are not due to heritable mutations of the SNPs, and that
  missingness is random and not e.g. due to heritable mutations of SNP
  flanking regions. Results based on this simulated data will provide an
  minimum estimate of the number of SNPs required, and an optimistic estimate
  of pedigree reconstruction performance.
}

\examples{
Geno_A <- SimGeno(Pedigree = Ped_griffin, nSnp=200, ParMis=c(0.1, 0.6),
                  MAF = 0.25, SnpError = 0.001)

Geno_B <- SimGeno(Pedigree = Ped_HSg5, nSnp = 100, ParMis = 0.2,
                 SnpError = c(0.01, 0.04, 0.1))

Geno_C <- SimGeno(Pedigree = Ped_griffin, nSnp=200, ParMis=0, CallRate=0.6,
                  SnpError = 0.05, ErrorFV=function(E) c(E/10, E/10, E))

# genotype matrix with duplicated samples:
Dups_grif <- data.frame(ID1 = c('i006_2001_M', 'i021_2002_M', 'i064_2004_F'))
Dups_grif$ID2 <- paste0(Dups_grif$ID1, '_2')
Err <- c(0.01, 0.04, 0.1)
Geno_act <- SimGeno(Ped_griffin, nSnp=500, ParMis=0, CallRate=1, SnpError=0)
Geno_sim <- MkGenoErrors(Geno_act, SnpError=Err, CallRate=0.99)
Geno_dups <- MkGenoErrors(Geno_act[Dups_grif$ID1, ], SnpError=Err,
                          CallRate=0.99)
rownames(Geno_dups) <- Dups_grif$ID2
Geno_sim <- rbind(Geno_sim, Geno_dups)

\dontrun{
# write simulated genotypes to a file, e.g. for use by PLINK:
GenoConvert(Geno_A, InFormat='seq', OutFormat='ped', OutFile = sim_genotypes)
}

}
\seealso{
The wrapper \code{\link{EstConf}} for repeated simulation and
  pedigree reconstruction; \code{\link{MkGenoErrors}} for fine control over
  the distribution of genotyping errors in simulated data;
  \code{\link{ErrToM}} for more information about genotyping error patterns.
}
\author{
Jisca Huisman, \email{jisca.huisman@gmail.com}
}
