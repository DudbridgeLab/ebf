#' Compute EBF bias for negative binomial tests
#' This is on hold as the integration over p seems problematical
#'
#' Computes the bias in the log posterior marginal likelihood for negative binomial tests.
#'
#' Pre-computed results are in the \code{ebf::nbinom.bias} object,
#' for sample sizes 1 to 100 with uniform prior.
#'
#' Full enumeration of the sample space can be very time-consuming.
#' Therefore the bias is estimated by Monte Carlo
#' simulations of observed and replicate data from the negative binomial distribution.
#'
#' @param size Vector of sample sizes for which to estimate bias.
#'
#' @param xmin Minimum probability of interval hypothesis.
#'
#' @param xmax Maximum probability of interval hypothesis.
#'
#' @param interval Number of probabilities to evaluate between \code{xmin} and \code{xmax}.
#'
#' @param shape Parameter of symmetrical Beta prior distribution.
#'
#' @param nsample Number of random samples taken to estimate bias.
#'
#' @param seed Random number seed
#'
#' @import stats
#'
#' @return A vector with the biases corresponding to the elements of
#' \code{size}.
#'
#' @template references
#'
#' @export

compute.nbinom.bias <- function(size=1:100, xmin=0, xmax=1, interval=100, shape=1,
                                complement=FALSE, nsample=1e6, seed=0) {

  # bias for each number of successes
  bias = NULL

  # loop through numbers of successes
  for(i in size) {

    # random number seed
    set.seed(seed)

    bias.i = NULL
    for(p in seq(xmin+1/interval, xmax-1/interval, length=interval)) {
      bias.p = NULL

      # observed data
      k.obs = rnbinom(nsample, i, p)
      if (complement == FALSE)
        b1 = pbeta(xmax, 2*i+shape, 2*k.obs+shape) -
        pbeta(xmin, 2*i+shape, 2*k.obs+shape)
      else
        b1 = pbeta(xmin, 2*i+shape, 2*k.obs+shape) +
        pbeta(xmax, 2*i+shape, 2*k.obs+shape, lower=F)

      # replicate data
      k.rep = rnbinom(nsample, i, p)
      if (complement == FALSE)
        b2 = pbeta(xmax, 2*i+shape, k.rep+k.obs+shape) -
        pbeta(xmin, 2*i+shape, k.rep+k.obs+shape)
      else
        b2 = pbeta(xmin, 2*i+shape, k.rep+k.obs+shape) +
        pbeta(xmax, 2*i+shape, k.rep+k.obs+shape, lower=F)

      # the bias for this probability
      w = which(b1 > 0 & b2 > 0)
      bias.p = lchoose(k.obs[w]+i-1, i-1) + lbeta(2*i+shape, 2*k.obs[w]+shape) + log(b1[w]) -
        lchoose(k.rep[w]+i-1, i-1) - lbeta(2*i+shape, k.rep[w]+k.obs[w]+shape) - log(b2[w])
      bias.i = c(bias.i, mean(bias.p))
    }
    bias = c(bias, mean(bias.i))
  }

  bias
}
