#' Compute EBF bias for binomial tests
#'
#' Computes the bias in the log posterior marginal likelihood for binomial tests.
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
  bias = NULL

  # loop through sample sizes
  for(i in n) {

    # find the integrated bias over the success probabilities
    bias = c(bias, integrate(function(p) {
      bias.i = 0
      # observed data
      for(k.obs in 0:i) {
        if (complement == FALSE)
          b1 = pbeta(xmax, 2*k.obs+shape, 2*(i-k.obs)+shape) -
            pbeta(xmin, 2*k.obs+shape, 2*(i-k.obs)+shape)
        else
          b1 = pbeta(xmin, 2*k.obs+shape, 2*(i-k.obs)+shape) +
            pbeta(xmax, 2*k.obs+shape, 2*(i-k.obs)+shape, lower=F)
        b1 = b1/dbinom(k.obs,i,0.5)
        # replicate data
        for(k.rep in 0:i) {
          if (complement == FALSE)
            b2 = pbeta(xmax, k.obs+k.rep+shape, 2*i-k.obs-k.rep+shape) -
              pbeta(xmin, k.obs+k.rep+shape, 2*i-k.obs-k.rep+shape)
          else
            b2 = pbeta(xmin, k.obs+k.rep+shape, 2*i-k.obs-k.rep+shape) +
              pbeta(xmax, k.obs+k.rep+shape, 2*i-k.obs-k.rep+shape, lower=F)
          b2 = b2/dbinom(k.rep,i,0.5)
          if (b1 > 0 & b2 > 0)
            bias.i = bias.i + dbinom(k.obs, i, p) * dbinom(k.rep, i, p) *
              (lchoose(i, k.obs) + lbeta(2*k.obs+shape, 2*(i-k.obs)+shape) + log(b1) -
                 lchoose(i, k.rep) - lbeta(k.obs+k.rep+shape, 2*i-k.obs-k.rep+shape) -
                 log(b2))
        }
      }

      # expectation over observed data
      bias.i^2 * dbeta(p, shape, shape)
    }, 0, 1)$value)
  }

  bias
}
