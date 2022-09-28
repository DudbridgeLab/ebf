#' @param h0 If a scalar, the value of a point null hypothesis.
#' If a vector with two elements, the lower and upper
#' bounds of the null hypothesis.
#' @param shrink If \code{TRUE}, uses the entire distribution of \code{x}
#' to calculate each individual EBF.
#'
#' If \code{shrink==TRUE} the calculation is of order \code{length(x)^2}
#' and may be very time-consuming.  To reduce the computation a subset of
#' \code{x} may be used to calculate the shrinkage EBFs.
#' By default a maximum of 1000 points is used, taken from the quantiles of the
#' simple EBFs.  This number can be varied with the \code{npoints} option.
#'
#' Further reduction can be achieved by calculating EBFs only for the elements
#' of \code{x} given by \code{index}.  The full distribution of \code{x} is used
#' in calculating those EBFs.
#'
#' @param npoints Number of points to use in calculate shrinkage EBFs.
#' Default is \code{min(length(x), 1000)}.
#'
#'
#' @return A data frame containing the following components:
#'   \itemize{
#'   \item{\code{index} Indices of elements in \code{x} for which EBFs are calculated.}
#'     \item{\code{ebf} Empirical Bayes factors for the elements of \code{x}.}
#'     \item{\code{ebf.units} Expressed as units of evidence, that is
#'      on the log scale with base 3.73.  One unit updates
#'      weaker belief to stronger belief, or weaker disbelief to weaker belief,
#'      according to a rational definition of weaker and stronger belief.
#'      See references for further details.}
#'     \item{\code{p} P-values.  This is empty if \code{h0} is an interval.}
#'     \item{\code{p.log10} Expressed in -log10 units.}
#'     \item{\code{ebf.shrink} Shrunken EBFs if specified by \code{shrink==TRUE}.}
#'     \item{\code{ebf.shrink.units} Expressed as units of evidence.}
#'   }
