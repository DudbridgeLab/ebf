#' Empirical Bayes factors for chi-squared tests
#'
#' Calculates empirical Bayes factors (EBFs) for chi-squared tests.
#'
#' @template allParams
#'
#' @param df Vector of degrees of freedom. Defaults to 1.
#'
#' @return A data frame containing the following components:
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
                      df=1,
                      index=NULL) {

  if (is.null(index)) index=1:length(x)

  ebf = exp((x-df)/2) / sqrt(2)^df
  ebf.units = log(ebf) / log((sqrt(3)+1)/(sqrt(3)-1))

  p = pchisq(x, df, lower=FALSE)
  p.log10 = -log(p)/log(10)

  result = data.frame(index =  index,
                      ebf = ebf[index],
                      ebf.units = ebf.units[index],
                      p = p[index],
                      p.log10 = p.log10[index])

  result
}
