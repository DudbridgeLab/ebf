#' Empirical Bayes factors for p-values
#'
#' Calculates empirical Bayes factors (EBFs) for stand-alone p-values.
#'
#' P-values are assumed to come from a Beta(1, b) distribution, with
#' a uniform prior on b from 1 to infinity.  Although this model is flexible,
#' there is no reason why it should hold exactly, and this EBF should
#' only be calculated when the p-value is the only available statistic.
#'
#' For small P-values, the EBF is approximately 1/10p.
#'
#' @template allParams
#' @template shrinkParams
#'
#' @import stats
#'
#' @template references
#'
#' @export

ebf.p <- function(x,
                  index=NULL,
                  shrink=FALSE,
                  points=NULL,
                  pi0=0) {

  if (is.null(index)) index=1:length(x)

  # EBFs
  ebf = ebf.p.simple(x)
  ebf.units = log(ebf) / log((sqrt(3)+1)/(sqrt(3)-1))

  # P-values
  p = x
  p.log10 = -log(p)/log(10)

  # shrinkage
  ebf.shrink = NULL
  ebf.shrink.units = NULL
  if (shrink) {

    # data points for estimating non-parametric distribution
    if (is.null(points)) points = 1:length(x)

    ebf.shrink = ebf.p.shrink(x, index, points, pi0)
    ebf.shrink.units = log(ebf.shrink) / log((sqrt(3)+1)/(sqrt(3)-1))
  }

  result = data.frame(index =  index,
                      ebf = ebf[index],
                      ebf.units = ebf.units[index])

  if (sum(!is.na(p[index])) > 0) result = data.frame(result,
                                                     p = p[index],
                                                     p.log10 = p.log10[index])

  if (shrink) result = data.frame(result,
                                  ebf.shrink = ebf.shrink[index],
                                  ebf.shrink.units = ebf.shrink.units[index])
  result
}
