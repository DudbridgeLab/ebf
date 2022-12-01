#' Compute EBF bias for Poisson tests
#'
#' Computes the bias in the log posterior marginal likelihood for Poisson tests.
#'
#' @param xmin Minimum probability of interval hypothesis.
#'
#' @param xmax Maximum probability of interval hypothesis.
#'
#' @param quant Minimum quantile of Poisson observations for which to compute bias.
#'
#' @param limit Upper limit of integration for Poisson rate.
#
#' @import stats
#'
#' @return The integrated bias over the rate parameter, which is 0.5
#'
#' @template references
#'
#' @export

compute.poisson.bias <- function(xmin=0, xmax=Inf, complement=FALSE, quant=1e-4, limit=100) {
  quant = min(quant, 1-quant)
  integrate(function(rate) {
    bias = 0
    for(i in min(qpois(quant, rate)):max(qpois(1-quant, rate))) {
      if (complement == FALSE)
        b1 = pgamma(xmax, 2*i+1, 2) - pgamma(xmin, 2*i+1, 2)
      else
        b1 = pgamma(xmin, 2*i+1, 2) + pgamma(xmax, 2*i+1, 2, lower=F)
      for(j in min(qpois(quant, rate)):max(qpois(1-quant, rate))) {
        if (complement == FALSE)
          b2 = pgamma(xmax, i+j+1, 2) - pgamma(xmin, i+j+1, 2)
        else
          b2 = pgamma(xmin, i+j+1, 2) + pgamma(xmin, i+j+1, 2, lower=F)
        if (b1 >0 & b2 > 0)
          bias = bias + dpois(i, rate) * dpois(j, rate) *
            ((j-i)*log(2) + lfactorial(2*i) + lfactorial(j) - lfactorial(i+j) - lfactorial(i) +
               log(b1) - log(b2))
      }
    }
    bias
  }, 0, limit)$value / limit

}
