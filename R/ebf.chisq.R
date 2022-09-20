#' Empirical Bayes factors for chi-squared tests
#'
#' Calculates empirical Bayes factors (EBFs) for chi-squared tests.
#'
#' @param x Vector of chi-squared statistics for which EBFs will be calculated.
#' @param df Vector of degrees of freedom. Defaults to 1.
#' @param index Vector of indices selecting a subset of \code{x}.
#'
#' @import stats
#'
#' @return An object of class "ebf" containing the following components:
#'   \itemize{
#'     \item{\code{ebf} Empirical Bayes factors for the elements of \code{x}.}
#'     \item{\code{ebf.units} \code{EBF} expressed as units of evidence.}
#'     \item{\code{p} P-values.}
#'     \item{\code{p.log10} -log10 of \code{p}.}
#'   }
#'
#' @author Frank Dudbridge
#'
#' @references
#' Dudbridge F (submitted) Units of evidence and expected Bayes factors for
#' objective reporting of statistical evidence.
#'
#' @export


ebf.chisq <- function(x,
                      df=NULL,
                      index=NULL) {

  if (is.null(df)) df = 1

  ebf = exp((x[index]-df)/2) / sqrt(2)^df
  ebf.units = log(ebf) / log((sqrt(3)+1)/(sqrt(3)-1))

  p = pchisq(x, df, lower=FALSE)
  p = p[index]
  p.log10 = -log(p)/log(10)

  result = list(index = index,
                ebf = ebf,
                ebf.units = ebf.units,
                p = p,
                p.log10 = p.log10)
  class(result) = ("ebf")
  result
}
