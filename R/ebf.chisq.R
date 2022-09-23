#' Empirical Bayes factors for chi-squared tests
#'
#' Calculates empirical Bayes factors (EBFs) for chi-squared tests.
#'
#' @param x Vector of chi-squared statistics for which EBFs will be calculated.
#' @param df Vector of degrees of freedom. Defaults to 1.
#'
#' @return An object of class "ebf" containing the following components:
#'   \itemize{
#'   \item{\code{index} Indices of elements in \code{x} for which EBFs are calculated.}
#'     \item{\code{ebf} Empirical Bayes factors for the elements of \code{x}.}
#'     \item{\code{ebf.units} Expressed as units of evidence.}
#'     \item{\code{p} P-values.  This is empty if \code{h0} is an interval.}
#'     \item{\code{p.log10} Expressed in -log10 units.}
#'   }
#'
#' @import stats
#'
#' @template references
#'
#' @export


ebf.chisq <- function(x,
                      df=NULL) {

  if (is.null(df)) df = 1

  ebf = exp((x-df)/2) / sqrt(2)^df
  ebf.units = log(ebf) / log((sqrt(3)+1)/(sqrt(3)-1))

  p = pchisq(x, df, lower=FALSE)
  p.log10 = -log(p)/log(10)

  result = list(index = index,
                ebf = ebf,
                ebf.units = ebf.units,
                p = p,
                p.log10 = p.log10)
  class(result) = ("ebf")
  result
}
