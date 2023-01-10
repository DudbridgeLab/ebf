#' Compute EBF bias for binomial tests
#'
#' Computes the bias in the log posterior marginal likelihood for binomial proportions.
#'
#' Pre-computed results are in the \code{ebf::binom.bias} object,
#' for sample sizes 1 to 100 with uniform prior.
#'
#' @param n Vector of sample sizes for which to compute bias.
#'
#' @param xmin Minimum probability of interval hypothesis.
#'
#' @param xmax Maximum probability of interval hypothesis.
#'
#' @param shape Parameter of symmetrical Beta prior distribution.
#'
#' @param complement If \code{TRUE}, computes the bias for the complement
#' of the interval hypothesis.
#'
#' @import stats
#'
#' @return A vector with the biases corresponding to the elements of
#' \code{n}.
#'
#' @template references
#'
#' @export

compute.binom.bias <- function(n=1:100, xmin=0, xmax=1, shape=1, complement=FALSE) {

  # bias for each sample size
  results = NULL

  # loop through sample sizes
  for(i in n) {

    bias = 0
    # observed data
    for(k.obs in 0:i) {
      if (complement == FALSE)
        b1 = pbeta(xmax, 2*k.obs+shape, 2*(i-k.obs)+shape) -
          pbeta(xmin, 2*k.obs+shape, 2*(i-k.obs)+shape)
      else
        b1 = pbeta(xmin, 2*k.obs+shape, 2*(i-k.obs)+shape) +
          pbeta(xmax, 2*k.obs+shape, 2*(i-k.obs)+shape, lower=F)

      # replicate data
      for(k.rep in 0:i) {
        if (complement == FALSE)
          b2 = pbeta(xmax, k.obs+k.rep+shape, 2*i-k.obs-k.rep+shape) -
            pbeta(xmin, k.obs+k.rep+shape, 2*i-k.obs-k.rep+shape)
        else
          b2 = pbeta(xmin, k.obs+k.rep+shape, 2*i-k.obs-k.rep+shape) +
            pbeta(xmax, k.obs+k.rep+shape, 2*i-k.obs-k.rep+shape, lower=F)

        if (b1 > 0 & b2 > 0) {
          bias = bias +
            # probability of k.obs and k.rep, marginal over p
            choose(i,k.obs) * choose(i, k.rep) *
            beta(k.obs+k.rep+shape, 2*i-k.obs-k.rep+shape) / beta(shape, shape) *
            # difference in log posterior marginal likelihoods
            (lchoose(i, k.obs) + lbeta(2*k.obs+shape, 2*(i-k.obs)+shape) + log(b1) -
               lbeta(k.obs+shape, i-k.obs+shape) -
               lchoose(i, k.obs) - lbeta(k.obs+k.rep+shape, 2*i-k.obs-k.rep+shape) -
               log(b2) +
               lbeta(k.rep+shape, i-k.rep+shape))
        }
      }
    }

   results = c(results, bias)
  }

  results
}
