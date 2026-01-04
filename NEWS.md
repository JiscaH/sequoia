# sequoia 3.2
- allow .vcf.gz as input for `GenoConvert`
- fix `which.min` applied to data.frame rows (no longer allowed) by adding `as.matrix`
- fix various issues affecting performance with hermaphrodites
- fix SeqOUT$DummyIDs output for hermaphrodite dummies
- hide documentation of internal functions


# sequoia 3.1.3
- add function `PlotPropAssigned`


# sequoia 3.1.2
- fix bug that disallowed sibship clustering among candidate parents in discrete 
2 generation offspring + candidate parents datasets.
- add function `CountOH`


# sequoia 3.1
- allow for age differences >100
- fix bug in `PedToNum` when called by `CalcPairLL` without pedigree
- `CalcPairLL`: when 'patmat' not specified, now always defaults to 'sex2', not
 only when focal=PO. 


# sequoia 3.0
- numerous updates to Fortran code to improve assignment accuracy, including several
  bug fixes related to monogamous mating systems
- fix inconsistency in ageprior matrix and its use for grandparents and avuncular pairs
- `PedCompare`: When dummy parent mismatched 1 or more dummy offspring than it had
  correct genotyped offspring, all genotyped + dummy offspring were set to 'nomatch'.
  Now match is kept, and only actual mismatched flagged as such. 
- `ComparePairs`: New argument 'Pairs_suffix', passed to `getRelM`. 
- `getRelM`: when only `Pairs` is provided, default (for unlisted pairs) changed
  from 'U' to 'X'
- fix bug in `SimGeno`, throwing an error when ParMis=c(0, not-zero)
- `CalcPairLL`: force to (almost) always calculate LL even when highly unlikely,
by fixing MaxMismatch to the total number of SNPs and T_filter to -999.0. 
- `CalcParentProbs`: new function to calculate assignment probabilities for any pedigree
- `getAssignCat` now also considers parents with only dummyfiable offspring as dummyfiable; 
 option '1sib' is dropped.
- `PedCompare`: when an incorrect half-sibling is added to a true singleton-sibship-with-grandparent,
  this now counts as 1 mismatch; this erroneously was 3 mismatches (entire sibship wrong).
- in `sequoia`, the default value for `CalcParentLLR` is changed from TRUE to FALSE.
- fix bug which caused no parents to be assigned to individuals produced by selfing
- add function `Err_Allele2Locus` 
 

# sequoia 2.11.5
- fix bug in Fortran code which causes some parents to not be assigned during `Module='ped'`
- fix 'Error in par(oldpar): invalid value specified for graphical parameter "plt"' when very small plotting window in newer versions of Rstudio. 

# sequoia 2.11.4
- `GenoConvert` fake .map output file: change chrom 0 (unmapped) to chrom 1, as SNPs on chromosome 0 get excluded by default by e.g. PLINK & GCTA
- `CalcRped`: drop dummy parents of 'half founders' from output (1 known + 1 unknown parent not supported by pkg kinship2)
- fix errors when using data.tables

# sequoia 2.11.3
- make assignment of grandparents to singletons a bit more conservative

# sequoia 2.11.2
- minor fixes to pass CRAN checks

# sequoia 2.11
- improved assignment rate when some birth years are unknown
- improved messages, with {cli} markup
- new format of runtime messages, using C script
- a few minor bugfixes in the Fortran code

# sequoia 2.9.1
- upgrade of `GenoConvert`: new vcf and genlight input, various bugs fixed, 
behaviour made more consistent and clearer through additional messages.
- retain SNP names in `GenoConvert` and `SnpStats`

# sequoia 2.9.0
- fix bug in ` MkGenoErrors` (used by `SimGeno`) causing about 3x too many hom|hom errors when `SnpError` is a single value: first beta-distributed per-SNP genotyping error rates 
  $E_l$ were generated, and then $(E_l/2)^2$ calculated. Now the single value is 
  by default first morphed into a length 3 vector (hom|hom, het|hom, hom|het), and three 
  beta distributions are generated. 
- adds log to `MkGenoErrors`
- `SimGeno` ParMis default changed from 0.4 to 0
- default genotyping pattern slightly changed to ensure the probability a homozyogote 
does not have a genotyping error is identical to a heterozygote (see `ErrToM`).
- Beta-version of `EstEr` (estimation of genotyping errors) removed due to inaccurate 
estimations and misuse. Will (probably) be re-implemented in a future version. 


# sequoia 2.8.3
- fixes bugs introduced since version 2.5, plus various other edits in source code to improve assignment rate
- fix bug in `CalcMaxMismatch`: OH with both parents counts as 2 mismatches (was 1)


# sequoia 2.7.3
- fix bug causing some negative parental LLRs, and possibly some non-assignments
- speed increase for lower call rates
- add OutFormat 'ped' to `GenoConvert`, and fix bug with OutFormat 'col'


# sequoia 2.7.2
- change `EstConf` example to nSim=1 to ensure runtime < 5 sec to pass CRAN check

# sequoia 2.7.1
- add `mtSame`: specify if individuals have the same or different mitochondrial haplotype
- improved parentage assignment performance when there are many genetically similar
candidate parents
- fixes CRAN issue 'cannot use Fortran's random number generator'


# sequoia 2.6.0
- add specification of assumed genotyping error rate via length 3 vector: hom|hom,
  het|hom, hom|het
- expand `CalcPairLL` helpfile 

# sequoia 2.5.6
- fixes CRAN pretest NOTES, including broken links in vignette

# sequoia 2.5.4
- add updated vignettes (main + age); accidentally included old versions in 2.5.3
- fixed bug in `CalcBYprobs`, which caused Year.last to be ignored

# sequoia 2.5.3
- fixes CRAN error 'DLL requires the use of native symbols'

# sequoia 2.5.1

### Bug fixes & minor changes
- fix error 'sibship number out of bounds'
- fix error when `LifeHistData$Sex` includes `NA`
- fix several minor bugs affecting rare cases


# sequoia 2.5.0

### New features & major changes
- improved performance when a large proportion of birth years is not exactly known
- optional column `Year.last` added to `LifeHistData` (last possible offspring birth year)
- New functions `GetAncestors` and `GetDescendants`

### Bug fixes & minor changes
- New parameter `MinAgeParent` for `MakeAgePrior()`
- New parameter `StrictGenoCheck`, `Strict` for `CheckGeno()` ; update msgs
- update `CheckLH`, now flexible column order in LifeHistData. 
- changed maxmismatch from `qntl = 0.999^(1/nrow(GenoM))` to `0.9999^(1/nrow(GenoM))` in 
all functions calling `CalcMaxMismatch`



# sequoia 2.4.2

### Bug fixes & minor changes
- checks up to 6 generations back when making assignment to avoid individual being its own ancestor (was 5)
- fix bug in `SnpStats()` when AF=0 or SNP is missing for all individuals. 


# sequoia 2.4.1

### New features & major changes
- R markdown file to create (or serve as a first draft of) a summary report of input and output, called via `sequoia_report()`
- `PedCompare()` parameter `minSibSize` (minimum sibship size for which the non-genotyped parent is considered 'dummyfiable') changed default value from `2sib` (2 genotyped siblings) to `1sib1GP`. This reflects the increased success of reconstructing grandoffspring-grandparent pairs in newest version, which make the output of `PedCompare(,minSibSize='2sib')` confusing. This also affects `EstConf()`. 
- `GetRelM()` now allows input of Pedigree AND Pairs. This a.o. allows `PlotRelPairs()` with both inferred pedigree plus `GetMaybeRel()` output.  



### Minor changes
- switch from library xlsx to openxlsx in `writeSeq()` (easier to install)
- adds purl=FALSE to vignette chunks & include `SeqOUT_HSg5` for quicker vignette compilation
- include additional griffin example data: `Geno_griffin`, `Conf_griffin` (output from `EstConf()`) and `MaybeRel_griffin` (output from `GetMaybeRel()`), as well as the script used to create these (`mk_griffin_data.R`)
- many examples are rewritten, for clarification & to speed up package check
- the plotting function of `SummarySeq()` was internal and is now exported (`PlotSeqSum()`)
- output from `sequoia()` now includes `args.AP`. 
- in `SummarySeq()`, if `Pedigree` is provided rather than `SeqList` and `SNPd=NULL`, all individuals are categorised as `Observed` (was: `Genotyped`).
- `PedPolish()` now has arguments to specify whether to drop extra columns (besides id-dam-sire) and whether to keep rows with non-unique or NA ids. 
- `SnpStats()` now includes HWE tests


### Bug fixes
- fixes several memory issues identified by CRAN
- fixes various bugs in `SimGeno()` for non-autosomal inheritance; this option is still experimental and non-autosomal SNPs are not supported by the pedigree reconstruction.
- fixes `cannot allocate vector of size ...` issue in `GetRelM()` with very large pedigrees, which affected `MakeAgePrior()` and thereby `sequoia()`. 
- `EstConf()` `$ConfProb` used the wrong denominator, namely the number of parents in the reference pedigree rather than in the inferred pedigree. 
- fixes bug resulting in LL(FA)=777 ('impossible') for pairs with all 4 parents unknown. Likely only affected `CalcPairLL()`, as LL(HS)=LL(GP)=LL(FA) for these pairs. Origin time unknown. 
- fixes bug in `GenoConvert()` when `Informat='single'`. 


# sequoia 2.3.5
- fixes bug: OH count was always zero when there was no co-parent


# sequoia 2.3.4
- fixes bug in `CalcPairLL()` HS likelihood when conditioning on pedigree was incorrect. No/minimal effect on pedigree reconstruction.
- fixes bug in `DuplicateCheck()` (always automatically called by `sequoia()`) that on very rare occasions caused R to crash


# sequoia 2.3.3
fixes minor bugs identified by CRAN valgrind and gcc-ASAN


# sequoia 2.3.2
minor edits to vignette to comply with CRAN precheck


# sequoia 2.3.1

### New features & major changes
- Hermaphrodites: dummy individuals with offspring as both dam and as sire now 
have prefix 'H'; closer links between the two 'clonal' sibship parts during pedigree reconstruction for improved performance
- Assignment of sibship grandparents moved to before check & assignment of additional parents; this proved to increase correct assignments without increasing incorrect assignments. 
- new function `CalcBYprobs()` to estimate the probability that individual i is born in year y. 

### Minor changes
- various edits to Fortran code improving general performance. 

### Bug fixes
- fixed inconsistent rounding in `EstConf()` output
- bug in `GenoConvert()` regarding InData vs InFile


# sequoia 2.3.0

### Major bug fixes
- now possible to run with assumed genotyping error rate `Err=0`
- further improvements when a large proportion of birth years are missing
- various bugs related to hermaphrodites

### Reduced computational time
- Sibship clustering etc. now done from oldest to youngest individual instead of in order of occurrence in genotype file, giving quicker convergence to final high likelihood pedigree.



# sequoia 2.2.2

### New features & major changes
- hermaphrodites: re-implemented and greatly improved sibship clustering. Specification of hermaphrodite vs diocious is now separate input parameter `Herm` instead of specified via `Complex`. New output list element `DummyClones`
- improved performance when a large proportion of birth years are missing
- new function `CalcRped()` to calculate pedigree relatedness. Uses package `kinship2`. 

### Bug fixes
Various smaller bugs have been fixed, some affecting assignment rate



# sequoia 2.1.0 

### New features & major changes
- parameter 'maxSibIter' (-9/-1/0/>0) deprecated and replaced by 'Module' (pre/dup/par/ped). 
- `sequoia()`: Option `FindMaybeRel` has been deprecated; call `GetMaybeRel()` directly instead.
- new function `CalcPairLL()`, returns likelihoods for each of the 7 considered relationships (PO, FS, HS, GP, FA, HA, U) for each specified pair of individuals
- new function `RelPlot()` for Colony-like visualisation of pairwise relationships (automatically called by `ComparePairs()`) 


### Bug fixes
- `sequoia()`: fixed error `object 'ErrM' not found` when re-using output from a previous sequoia run. Circumvent this bug in version 2.0.7 by fooling the program to think it's output from an older version: `names(ParOUT$Specs)[match("MaxMismatchOH", names(ParOUT$Specs))] <- "foo"`. 
- `sequoia()`: fixed bug causing genotyped parents to not always be monogamous when `Complx='mono'`.
- various functions: fixed error when dummy prefixes have different number of characters (`Error in data.frame(id = c(dID[s(nd[1]), 1], dID[s(nd[2]), 2]), VtoM(TMP$dumparrf,  : arguments imply differing number of rows)` )
- `GetMaybeRel()`: fixed error `(subscript) logical subscript too long` when input pedigree contains dummies
- `GetMaybeRel()`: fixed error causing likely GP pairs not to be included output
- `PedCompare()`: fixed id.dam.cat and id.sire.cat being 'NANA' instead of XD, XG or XX when Symmetrical=TRUE
- `PlotAgePrior()`: Avoid using `grDevices::hcl.colors()` in R versions <3.6, where this function is not yet available 
- `ComparePairs()`: fixed bug when `Pairs2` but not `Ped2` is specified
- `MakeAgePrior()`: fixed bug when there are no FS pairs in the input pedigree
- `MakeAgePrior()`: `MaxAgeParent` was ignored when a pedigree with overlapping generations was supplied


### Minor changes
- most examples are now set to `\donttest` instead of `\dontrun`, so that they can be run using `example()`. Note that some can be quite time consuming, especially `EstConf()`.
- new function `PlotPedComp()` to visualise `PedCompare()` output
- `SimGeno()`: deprecated input parameters (since v 1.3.1) are dropped completely
- `getAssignCat()` no longer drops additional columns from input pedigree
- Speed increase for `CalcOHLLR(, CalcLLR = FALSE)`
- More thorough input checks, which are more consistent across different functions
- `PedCompare()` output element `DummyMatch` now also include output class of matched individual's parents & offspring
- Duplicate check in `sequoia()` now only returns pairs for which LL_duplicate - max(LL_{not duplicate}) > T_filter; when call rates are low this may be a substantially shorter list than in previous versions, where all pairs with fewer than MaxMismatchDUP differences were listed. 
- `ComparePairs()` can now be called for a single pedigree, as well as to compare two pedigrees
- If the plotting window in Rstudio is too small and `Plot=TRUE`, all functions will print a message and return results as usual, instead of throwing an error and not returning results.
- `MakeAgePrior()`: more consistent implementation; is now called by `sequoia()` with only  lifehistory data of genotyped individuals. 
- `EstConf()` now also returns the full `Counts` table from `PedCompare()`; `$RunParams` now holds the evaluated input paramters, instead of e.g. `V[i]` when called from inside a loop. 


# sequoia 2.0.7 

### New features
- The genotyping error matrix (probability of observed genotype conditional on actual) is now fully customisable in all relevant functions, see help file of new function `ErrToM`. The default has changed very slightly from version 1.3. 
- new function `CalcOHLLR()` to calculate Mendelian errors and parental log-likelihood ratios for any pedigree
- new function `getAssignable()` to flag genotyped and 'dummifiable' individuals in any pedigree
- new function `ComparePairs()` to compare pairwise relationships between 2 pedigrees; replaces now-deprecated `DyadCompare`. 
- function `PedPolish()` is now user available.

### Major changes
- Deprecated option `MaxMismatch` of function `sequoia`, now calculated internally by new function `CalcMaxMismatch` based on number of SNPs, presumed genotyping error rate, and minor allele frequencies
- function `EstConf()` now also estimates confidence for parent-pairs; output has changed considerably. 
- rewrote function `PedCompare()` to increase clarity of code for easier maintenance; changed output format somewhat. 
- Added a vignette about the ageprior, and rewrote sections of the main vignette to incorporate new functions
- In the Fortran part, re-implemented how (candidate) (grand)parent-pairs are filtered and assigned 


### Bug fixes
- `GenoConvert()` skipped first individual when reading .raw file. Circumvent this bug in earlier versions by using option `header=FALSE` (then header row is removed only once...)
- `ConfProb()` expected input parameter `nSim` to be strictly integer, now relaxed to any value convertible to a whole number
- fixed `ERROR! ***Invalid ParProb!***` triggered when some SNPs are monomorphic
- fixed `SEGFAULT` triggered when some SNPs have very high missingness (>80%); possibly sibship size out of bounds
- fixed `Error arguments imply differing number of rows` when there are dummy parents of 1 sex only
- fixed various mostly minor bugs in the Fortran code
- fixed bugs regarding 'link time optimisation'



### Minor changes
- `LifeHistData` may have 2 additional columns, with minimum and maximum possible birth year
- second example pedigree (`Ped_griffin`) to illustrate overlapping generations, used in age vignette 
- `SummarySeq()`: added a pedigree summary table identical to a subset of the table returned by R package `pedantics`' `pedStatSummary`; that package has been archived on CRAN. Added option `Panels` to only plot (a) specific panel(s).  


# Sequoia 1.3.3

### Bug fixes
- fixes bug that caused R to crash (Fortran array indexing out-of-bounds)


# Sequoia 1.3.1

### New features
- several functions have become user-visible: `CheckGeno()`, `MkGenoErrors()`, `GetMaybeRel()`, `GetRelCat()`
- plotting functions added: `PlotAgePrior()` and `SummarySeq()`
- function `SimGeno()` input parameters have changed, old ones will be deprecated


### Minor changes
- extended vignette with function overview & FAQ
- numerous edits to fortran source code to better handle certain (rarer) types of relatives


### Bug fixes
- various bug fixes in fortran source code


# sequoia 1.1.1

### Major changes
- possibly. 


# sequoia 1.0.0

### New features
- added functions `EstConf`, `SnpStats`

### Major changes
- considerable changes in Fortran code


# sequoia 0.9.2

### New features
- added functionality for hermaphrodites (in silico cloned into male + female)


# sequoia 0.7.2
First version on CRAN!
