#' @title Plot Age Priors
#'
#' @description Visualise the age-difference based prior probability ratios as a
#'   heatmap.
#'
#' @param AP  matrix with age priors (\eqn{P(A|R)/P(A)}) with age differences in
#'   rows and relationships in columns; by default M: maternal parent (mother),
#'   P: paternal parent (father), FS: full siblings, MS: maternal siblings (full
#'   + half), PS: paternal siblings.
#' @param legend if \code{TRUE}, a new plotting window is started and
#'   \code{\link{layout}} is used to plot a legend next to the main plot. Set to
#'   \code{FALSE} if you want to add it as panel to an existing plot (e.g. with
#'   \code{par(mfcol=c(2,2))}).
#'
#' @return A heatmap.
#'
#' @seealso \code{\link{MakeAgePrior}}, \code{\link{SummarySeq}}.
#'
#' @importFrom graphics layout par image axis mtext abline
#' @importFrom grDevices hcl
#'
#' @examples
#' PlotAgePrior(SeqOUT_griffin$AgePriors)
#' PlotAgePrior(SeqOUT_griffin$AgePriorExtra)
#'
#' @export

PlotAgePrior <- function(AP = NULL, legend=TRUE)
{
  if (is.data.frame(AP))  AP <- as.matrix(AP)
  if (!is.matrix(AP))  stop("AP must be a matrix")
  if (any(AP < 0 | AP > 1000) | any(!is.double(AP)))
    stop("AP must be a numeric matrix with values >0 & < 1000")
  RR <- colnames(AP)
  if (ncol(AP)==15 && all(c("M", "P", "FS", "MS", "PS", "MGF", "PGM", "MFA", "PPA") %in% RR)) {  # re-arrange for clarity
    RR <- c("M", "P", "FS", "MS", "PS",
        "MGM", "MGF", "PGM", "PGF",
        "MFA", "PFA", "MMA", "MPA", "PMA", "PPA")
  }

  # trim off rows with only zero's. First row: A=0
  MaxA <- min(max(which(apply(AP, 1, function(x) any(x>0)))) +1, nrow(AP))
  if (!is.null(rownames(AP))) {
    AA <- as.numeric(rownames(AP))[1:MaxA]
  } else {
    AA <- 1:MaxA
  }

	stretch <- function(V, x=10) {
  	a <- list()
  	for (i in 1:(length(V)-1)) {
  	  a[[i]] <- seq(V[i], V[i+1], length.out=x+1)[1:x]
  	}
	  return(c(unlist(a), V[length(V)]))
	}
	flip <- function(M, dim=2)  apply(M, dim, function(x) rev(x))

	#~~~  heatmap  ~~~~~~~
	brks <- c(1:5, 10,20,50)
	brks <- c(0, 1/1001, rev(1/brks[-1]), brks) - 1e-5
	brks.f <- stretch(brks, x=10)
	mids <- (brks.f[-length(brks.f)] + brks.f[-1])/2
	if (exists("hcl.colors")) {
	  cols <- c(1, grDevices::hcl.colors(89, "Light Grays"),
	            grDevices::hcl.colors(70, "Greens 3", rev=TRUE))
	} else {
	  # colour specification obtained via colorspace::hcl_palettes()
	  seqp <- function(from, to, length.out, pwr) {
	    i <- seq.int(1, 0, length.out = length.out)
	    from + (to - from) * (1 - i^pwr)
	  }
	  cols <- c(1, hcl(h=0, c=0, l=seqp(30, 90, len=89, pwr=1.5)),
	            rev(hcl(h=135,
	                c=seq(100, 2, len=70),
	                l=seqp(25, 98, len=70, pwr=1.5))))
	}

	oldpar <- par(no.readonly = TRUE)
	oldpar <- oldpar[!names(oldpar) %in% c("pin","fig",'plt')]   # current plot dimensions, not setable. bug?

	if (legend) {
	  ly <- tryCatch( layout(matrix(c(1,2), nrow=1), widths=c(.8, .2)) ,
	                  error = function(e) {
	                    message("PlotAgePrior: Plotting area too small for legend")
	                    return(NA)
	                  } )
	  if (is.na(ly)) {
	    legend <- FALSE
	  } else {
	    par(mai=c(.9,.9,.2,.1))
	  }
	}
	img <- tryCatch(
	  {
	    suppressWarnings(image(x=t(AP[1:MaxA, RR]), y=AA, xaxt="n",
	                           yaxt=ifelse(length(AA)<5, "n", "s"), las=1,
	                           breaks = brks.f, col = cols,
	                           ylab="Age difference (A)", cex.lab=1.1))
	  },
	  error = function(e) {
	    message("PlotAgePrior: Plotting area too small for heatmap")
	    return(NA)
	  })
	if (!is.null(img)) {
	  par(oldpar)
	  return()
	}

  if (length(AA)<5)  axis(side=2, at=AA, labels=AA, las=1, cex.lab=1.1)
	axis(side=1, at=seq(0,1,along.with=RR), line=-0.5, labels=RR,
	     cex.axis=ifelse(length(RR)==5, 1.1, 0.8),
	     las = ifelse(length(RR)==5, 1, 2),
	     tick=FALSE)
	if (any(AA<0))  abline(h=0, col=2)  #abline(h=-0.4, col=2)
	mtext("Relationship (R)", side=1, line=2, cex=1.1)

	if (legend) {
  	par(mai=c(1.1,.3,0.5,.7))  # legend
	  E <- tryCatch(suppressWarnings(image(t(as.matrix(mids)), axes=FALSE,
	                                       frame.plot=TRUE, breaks = brks.f, col = cols)),
	                error = function(e) {
	                  message("PlotAgePrior: Plotting area too small for legend")
	                  return(NA)
	                } )
  	if (is.null(E)) {  # sometimes doesn't fit in plotting window
    	axis(side=4, at=seq(0, 1,length.out=length(brks)), labels=round(brks, 3), las=1, cex.axis=0.8)
    	axis(side=4, at=seq(0, 1,length.out=length(brks)), labels=FALSE, tck=0.2) # tcl=0.5)
      axis(side=2, at=seq(0, 1,length.out=length(brks)), labels=FALSE, tck=0.2)
    	mtext(("P(A|R)/P(A)"), side=3, line=0.5, cex=1)
  	} else {
  	  par(mai=rep(0,4))
      graphics::plot.new()  # so that can continue normally
  	}
	  par(oldpar)  # restore old par settings  (if !legend, par not changed)
	}
}
