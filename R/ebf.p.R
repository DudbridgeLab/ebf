#' Empirical Bayes factors for p-values
#'
#' Calculates empirical Bayes factors (EBFs) for standalone p-values.
#'
#' P-values are assumed to come from a Beta(1, b) distribution, with
#' a uniform prior on b from 1 to infinity.  Although this model is flexible,
#' there is no reason why it should hold exactly, and this EBF should
#' only be calculated when the p-value is the only available statistic.
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
                  npoints=1000,
                  nsupport=20,
                  tol=1e-5,
                  nboot=0,
                  seed=0) {

  set.seed(seed)
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
    if (length(x) < npoints) points = 1:length(x)
    else points = sample(1:length(x), npoints)

    # number of support points in non-parametric distribution
    nsupport = min(length(x), nsupport)
    ebf.shrink = ebf.p.npml(x, index, points, nsupport, tol, nboot)
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
