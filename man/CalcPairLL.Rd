% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CalcPairLL.R
\name{CalcPairLL}
\alias{CalcPairLL}
\title{Calculate Likelihoods for Alternative Relationships}
\usage{
CalcPairLL(
  Pairs = NULL,
  GenoM = NULL,
  Pedigree = NULL,
  LifeHistData = NULL,
  AgePrior = TRUE,
  SeqList = NULL,
  Module = "ped",
  Complex = "full",
  Herm = "no",
  InclDup = FALSE,
  Err = 1e-04,
  ErrFlavour = "version2.9",
  Tassign = 0.5,
  Tfilter = -2,
  quiet = FALSE,
  Plot = TRUE
)
}
\arguments{
\item{Pairs}{dataframe with columns \code{ID1} and \code{ID2}, and
 optionally
 \describe{
\item{Sex1}{Sex of ID1, 1=female, 2=male, 3=unknown, or \code{NA} to take
  from \code{LifeHistData}. The sex of individuals occurring as parent in
  \code{Pedigree} cannot be altered.}
\item{Sex2}{Sex of ID2}
\item{AgeDif}{Age difference in whole time units, BirthYear1 - BirthYear2
  (i.e. positive if ID2 is born before ID1). If \code{NA}, calculated from
  \code{LifeHistData}. Use '999' to explicitly specify 'unknown'.}
\item{focal}{relationship character abbreviation; PO, FS, HS, GP or U. See
Details for its effect and explanation of abbreviations. Default: U.}
\item{patmat}{1=maternal relatives, 2=paternal relatives. Only relevant for
  HS & GP, for which it defaults to Sex1, or 1 if Sex1=3, but is currently
  only predictably implemented for pairs of two genotyped individuals.
  Always equal to Sex2 for PO pairs when Sex2 is known.}
\item{dropPar1}{Drop the parents of \code{ID1} before calculating the pair
likelihood, rather than conditioning on them; choose from 'none', 'dam',
'sire', or 'both'. See example. If e.g. the pair shares a common mother,
'none' and 'sire' will condition on this shared mother and not calculate the
likelihood that they are maternal siblings, while dropPar1='dam' or 'both'
will calculate that likelihood, and the other likelihoods as if the mother
of ID1 were unknown.}
\item{dropPar2}{as \code{dropPar1}, for \code{ID2}}
}}

\item{GenoM}{numeric matrix with genotype data: One row per individual,
one column per SNP, coded as 0, 1, 2, missing values as a negative number
or NA. You can reformat data with \code{\link{GenoConvert}}, or use other
packages to get it into a genlight object and then use \code{as.matrix}.}

\item{Pedigree}{dataframe with columns id-dam-sire; likelihoods will be
calculated conditional on the pedigree. May include non-genotyped
individuals, which will be treated as dummy individuals.}

\item{LifeHistData}{data.frame with up to 6 columns:
 \describe{
 \item{ID}{max. 30 characters long}
 \item{Sex}{1 = female, 2 = male, 3 = unknown, 4 = hermaphrodite,
           other numbers or NA = unknown}
\item{BirthYear }{birth or hatching year, integer, with missing values as NA
  or any negative number.}
\item{BY.min}{minimum birth year, only used if BirthYear is missing}
\item{BY.max}{maximum birth year, only used if BirthYear is missing}
\item{Year.last}{Last year in which individual could have had offspring. Can
  e.g. in mammals be the year before death for females, and year after death
  for males. } }
"Birth year" may be in any arbitrary discrete time unit relevant to the
species (day, month, decade), as long as parents are never born in the same
time unit as their offspring, and only integers are used. Individuals do not
need to be in the same order as in `GenoM', nor do all genotyped individuals
need to be included.}

\item{AgePrior}{logical (\code{TRUE/FALSE}) whether to estimate the ageprior
from \code{Pedigree} and \code{LifeHistData}, or a matrix as generated by
\code{\link{MakeAgePrior}} and included in the \code{\link{sequoia}}
output. The \code{AgePrior} affects which relationships are considered
possible: only those where \eqn{P(A|R) / P(A) > 0}.  When \code{TRUE},
\code{\link{MakeAgePrior}} is called using its default values. When
\code{FALSE}, all relationships are considered possible for all age
differences, except that parent-offspring pairs cannot have age difference
zero, and grand-parental pairs have an age difference of at least two.}

\item{SeqList}{list with output from \code{\link{sequoia}}. If input
parameter \code{Pedigree=NULL}, \code{SeqList$Pedigree} will be used if
present, and \code{SeqList$PedigreePar} otherwise. If \code{SeqList$Specs}
is present, input parameters with the same name as its items are ignored.
The list elements 'LifeHist', 'AgePriors', and 'ErrM' are also used if
present, and override the corresponding input parameters.}

\item{Module}{if \code{ped} (full pedigree), turn any non-genotyped parents
in \code{Pedigree} into dummies to condition upon. If \code{par}, ignore
any non-genotyped parents.}

\item{Complex}{Breeding system complexity. Either "full" (default), "simp"
(simplified, no explicit consideration of inbred relationships), "mono"
(monogamous).}

\item{Herm}{Hermaphrodites, either "no", "A" (distinguish between dam and
sire role, default if at least 1 individual with sex=4), or "B" (no
distinction between dam and sire role). Both of the latter deal with
selfing.}

\item{InclDup}{logical, include the likelihood for the two samples to be
duplicates (originating from the same individual) in the output?}

\item{Err}{estimated genotyping error rate, as a single number, or a length 3
vector with P(hom|hom), P(het|hom), P(hom|het), or a 3x3 matrix. See
details below. The error rate is presumed constant across SNPs, and
missingness is presumed random with respect to actual genotype. Using
\code{Err} >5\% is not recommended, and \code{Err} >10\% strongly
discouraged.}

\item{ErrFlavour}{function that takes \code{Err} (single number) as input,
and returns a length 3 vector or 3x3 matrix, or choose from inbuilt options
'version2.9', 'version2.0', 'version1.3', or 'version1.1', referring to the
sequoia version in which they were the default. Ignored if \code{Err} is a
vector or matrix. See \code{\link{ErrToM}} for details.}

\item{Tassign}{minimum LLR required for acceptance of proposed relationship,
relative to next most likely relationship. Higher values result in more
conservative assignments. Must be zero or positive.}

\item{Tfilter}{threshold log10-likelihood ratio (LLR) between a proposed
relationship versus unrelated, to select candidate relatives. Typically a
negative value, related to the fact that unconditional likelihoods are
calculated during the filtering steps. More negative values may decrease
non-assignment, but will increase computational time.}

\item{quiet}{logical, suppress messages}

\item{Plot}{logical, display scatter plots by \code{\link{PlotPairLL}}.}
}
\value{
The \code{Pairs} dataframe including all optional columns listed
  above, plus the additional columns:
 \item{xx}{Log10-Likelihood of this pair having relationship xx, with xx
   being the relationship abbreviations listed below.}
 \item{TopRel}{Abbreviation of most likely relationship}
 \item{LLR}{Log10-Likelihood ratio between most-likely and second most likely
   relationships. Other LLRs, e.g. between most-likely and unrelated, can
   easily be computed.}

\strong{Relationship abbreviations:}
  \item{PO}{Parent - offspring}
  \item{FS}{Full siblings}
  \item{HS}{Half siblings}
  \item{GP}{Grandparent}
  \item{FA}{Full avuncular}
  \item{HA}{Half avuncular and other 3rd degree relationships}
  \item{U}{Unrelated}
  \item{2nd}{Unclear which type of 2nd degree relatives
    (HS, GP, or FA)}
  \item{??}{Unclear which type of 1st, 2nd or 3rd degree
    relatives}

\strong{Likelihood special codes:}
  \item{222}{Maybe (via) other parent (e.g. focal="GP", but as likely to be
    maternal as paternal grandparent, and therefore not assignable)}
  \item{333}{Excluded from comparison (shouldn't occur)}
  \item{444}{Not implemented (e.g. would create an odd double/triple
    relationship in combination with the provided pedigree)}
  \item{777}{Impossible (e.g. a male (Sex2=2) cannot be mother (patmat=1))}
  \item{888}{Already assigned in the provided pedigree (see \code{dropPar}
    arguments)}
  \item{999}{\code{NA}. If all values for the pair are \code{999}, one or both
    individuals are not genotyped and not dummyfiable.}
}
\description{
For each specified pair of individuals, calculate the
  log10-likelihoods of being PO, FS, HS, GP, FA, HA, U (see Details).
  Individuals must be genotyped or have at least one genotyped offspring.

  \strong{NOTE} values \eqn{>0} are various \code{NA} types, see 'Likelihood
  special codes' in 'Value' section below.
}
\details{
The same pair may be included multiple times, e.g. with different
  sex, age difference, or focal relationship, to explore their effect on the
  likelihoods. Likelihoods are only calculated for relationships that are
  possible given the age difference, e.g. PO (parent-offspring) is not
  calculated for pairs with an age difference of 0.

  Non-genotyped individuals can be included if they have at least one
  genotyped offspring and can be turned into a dummy (see
  \code{\link{getAssignCat}}); to establish this a pedigree must be provided.

  \strong{Warning 1}: There is no check whether the input pedigree is genetically
  sensible, it is simply conditioned upon. Checking whether a pedigree is
  compatible with the SNP data can be done with \code{\link{CalcOHLLR}}.

  \strong{Warning 2}: Conditioning on a \code{Pedigree} can make computation
  orders of magnitude slower.
}
\section{Why does it say 444 (no can do) for all relationships?}{

 This happens when the pair does not pass the initial check which prevents
 impossible configurations in combination with \code{Pedigree}. Specifically,
 it happens when either or both individuals are a parent in the pedigree, but
 \itemize{
   \item the sex in \code{Pairs} is not consistent with that
   \item \code{Pairs} changes the age differences; it is too complex to check
   whether or not this still makes all pedigree links valid. Only setting the
   age difference to 'unknown' via \code{Pairs} is possible.}
}

\section{Why does it say 777 (impossible)?}{

  This function uses the same machinery as \code{sequoia}, which will to save
  time not calculate the likelihood when it is quickly obvious that the pair
  cannot be related in the specified manner.

  For PO (putative parent-offspring pairs) this is the case when:
  \itemize{
   \item the sex of the candidate parent, via \code{Pairs$Sex2} or
    \code{LifeHistData}, does not match \code{Pairs$patmat}, which defaults
    to 1 (maternal relatives, i.e. dam)
  \item a dam is already assigned via \code{Pedigree} and \code{Pairs$dropPar1
    ='none'}, and \code{Pairs$patmat = 1}
  \item the age difference is zero or otherwise impossible according to the
    age prior. It is either calculated from \code{LifeHistData} or specified
    via \code{Pairs$AgeDif}. The \code{AgePrior} can be specified directly,
    be taken from \code{SeqList}, or calculated  automatically by
    \code{\link{MakeAgePrior}} when both \code{Pedigree} and
    \code{LifeHistData} are provided. }
}

\section{Double relationships & focal relationship}{

  Especially when Complex='full', not only the seven relationship
  alternatives listed above are considered, but a whole range of possible
  double and even triple relationships. For example, mother A and offspring B
  (PO) may also be paternal half-siblings (HS, A and A's mother mated with
  same male), grandmother and grand-offspring (GP, B's father is A's son), or
  paternal aunt (B's father is a full or half sib of A).

  The likelihood reported as 'LL_PO' is the most-likely one of the possible
  alternatives, among those that are not impossible due to age differences or
  due to the pedigree (as reconstructed up to that point). Whether e.g. the
  likelihood to be both PO & HS is counted as PO or as HS, depends on the
  situation and is determined by the variable 'focal': During parentage
  assignment, it is counted as PO but not HS, while during sibship
  clustering, it is counted as HS but not PO -- not omitting from the
  alternative relationship would result in a deadlock.
}

\section{Dummy individuals}{

  For historical reasons, the relationships between a dummy ID1 and ID2 are
  reported *between the sibship and ID2*. So,
  \describe{
 \item{PO}{ID2 replaces dummy ID1; or merge dummy ID2 with dummy ID1}
 \item{FS, HS}{ID1 parent of ID2}
 \item{GP}{ID2 parent of ID1}
 \item{FA,HA}{ID2 FS resp. HS of ID1} }
 If ID1 is genotyped and ID2 is a dummy, the relationships are as when ID2 is
 genotyped.

 If you wish to retrieve likelihoods for a different set of relationships,
 please contact me at \email{jisca.huisman@gmail.com} .
}

\examples{
# Likelihoods depend on the presumed genotyping error rate:
CalcPairLL(Pairs = data.frame(ID1='i042_2003_F', ID2='i015_2001_F'),
           GenoM = Geno_griffin, Err = 1e-7, Plot=FALSE)
CalcPairLL(Pairs = data.frame(ID1='i042_2003_F', ID2='i015_2001_F'),
           GenoM = Geno_griffin, Err = 1e-3, Plot=FALSE)


## likelihoods underlying parent LLR in pedigree:
# Example: dams for bottom 3 individuals
tail(SeqOUT_griffin$PedigreePar, n=3)
# set up dataframe with these pairs. LLRdam & LLRsire ignore any co-parent
Pairs_d <- data.frame(ID1 = SeqOUT_griffin$PedigreePar$id[140:142],
                      ID2 = SeqOUT_griffin$PedigreePar$dam[140:142],
                      focal = "PO",
                      dropPar1 = 'both')

# Calculate LL's, conditional on the rest of the pedigree + age differences
CalcPairLL(Pairs_d, GenoM = Geno_griffin, Err = 1e-04,
           LifeHistData = LH_griffin, Pedigree = SeqOUT_griffin$PedigreePar)

# LLR changes when ignoring age and/or pedigree, as different relationships
# become (im)possible
CalcPairLL(Pairs_d, GenoM = Geno_griffin, Err = 1e-04)

# LLRpair is calculated conditional on co-parent, and min. of dam & sire LLR
Pairs_d$dropPar1 <- 'dam'
Pairs_s <- data.frame(ID1 = SeqOUT_griffin$PedigreePar$id[141:142],
                      ID2 = SeqOUT_griffin$PedigreePar$sire[141:142],
                      focal = "PO",
                      dropPar1 = 'sire')
CalcPairLL(rbind(Pairs_d, Pairs_s), GenoM = Geno_griffin, Err = 1e-04,
           LifeHistData = LH_griffin, Pedigree = SeqOUT_griffin$PedigreePar)


## likelihoods underlying LLR in getMaybeRel output:
MaybeRel_griffin$MaybePar[1:5, ]
FivePairs <- MaybeRel_griffin$MaybePar[1:5, c("ID1", "ID2", "Sex1", "Sex2")]
PairLL <- CalcPairLL(Pairs = rbind( cbind(FivePairs, focal = "PO"),
                                    cbind(FivePairs, focal = "HS"),
                                    cbind(FivePairs, focal = "GP")),
                     GenoM = Geno_griffin, Plot=FALSE)
PairLL[PairLL$ID1=="i121_2007_M", ]
# LL(FS)==222 : HSHA, HSGP, and/or FAHA more likely than FS
# LL(GP) higher when focal=HS: GP via 'other' parent also considered
# LL(FA) higher when focal=PO: FAHA, or FS of 'other' parent

}
\seealso{
\code{\link{PlotPairLL}} to plot alternative relationship pairs from
  the output; \code{\link{LLtoProb}} to transform likelihoods to
  probabilities; \code{\link{CalcParentProbs}} which uses this function to
  calculate parental probabilities; \code{\link{GetRelM}} to find all
  pairwise relatives according to the pedigree; \code{\link{GetMaybeRel}} to
  get likely relative pairs based on the genetic data.
}
