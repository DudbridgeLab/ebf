#' @param x Vector of statistics for which EBFs will be calculated.
#' @param se Vector of standard errors associated with \code{x}.  Default is 1.
#' @param h0 A vector with two elements, giving the lower and upper
#' bounds of the null hypothesis.
#' @param h1 A vector with two elements, giving the lower and upper bounds of
#' the alternative hypothesis.  If unspecified, the alternative
#' hypothesis is the complement of \code{h0}.
#' @param shrink If \code{TRUE}, uses the entire distribution of \code{x}
#' to calculate each individual EBF.
#' @param npoints Number of points to use in calculate shrinkage EBFs.
#' Default is \code{min(length(x), 1000)}.
#' @param index Vector of indices of elements of \code{x} for which to calculate shrinkage EBFs.
#'
#' @return An object of class "ebf" containing the following components:
#'   \itemize{
#'   \item{\code{index} Indices of elements in \code{x} for which EBFs are calculated.}
#'     \item{\code{ebf} Empirical Bayes factors for the elements of \code{x}.}
#'     \item{\code{ebf.units} Expressed as units of evidence.}
#'     \item{\code{ebf.shrink} Shrunken EBFs if specified by \code{shrink==TRUE}.}
#'     \item{\code{ebf.shrink.units} Expressed as units of evidence.}
#'     \item{\code{p} P-values.  This is empty if \code{h0} is an interval.}
#'     \item{\code{p.log10} Expressed in -log10 units.}
#'   }
#
