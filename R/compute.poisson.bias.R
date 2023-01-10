#' Compute EBF bias for Poisson tests
#'
#' Computes the bias in the log posterior marginal likelihood for Poisson rates.
#'
#' @param xmin Minimum probability of interval hypothesis.
#'
#' @param xmax Maximum probability of interval hypothesis.
#'
#' @param shape Shape parameter of prior Gamma distribution.
#'
#' @param rate Rate parameter of prior Gamma distribution.
#'
#' @param quant Minimum quantile of Poisson observations for which to compute bias.
#'
#' @param limit Upper limit of observed data.
#
#' @import stats
#'
#' @return The bias, which is 0.5 when \code{xmax==Inf} and 0 otherwise.
#'
#' @template references
#'
#' @export

compute.poisson.bias <- function(xmin=0, xmax=Inf, shape=1, rate=0,
                                 quant=1e-4, limit=100) {

  quant = min(quant, 1-quant)

  bias = 0
  for(i in 0:limit) {
    b1 = pgamma(xmax, 2*i+shape, 2+rate) - pgamma(xmin, 2*i+shape, 2+rate)
    for(j in qpois(quant, qgamma(quant, i+1, 1)):qpois(1-quant, qgamma(1-quant, i+1, 1))) {
      b2 = pgamma(xmax, i+j+shape, 2+rate) - pgamma(xmin, i+j+shape, 2+rate)
      if (b1 >0 & b2 > 0)
        bias = bias +
          # probability of j and i, marginal over lambda
          exp(lfactorial(i+j) - (i+j+1)*log(2) - lfactorial(i) - lfactorial(j)) /
          limit *
          # difference in log posterior marginal likelihoods
          ((j-i)*log(2+rate) + lgamma(2*i+shape) + lgamma(j+shape) -
             lgamma(i+j+shape) - lgamma(i+shape) +
             log(b1) - log(b2))
    }
  }

  bias
}
