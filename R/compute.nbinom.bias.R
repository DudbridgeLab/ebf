#' Compute EBF bias for negative binomial tests
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

compute.nbinom.bias <- function(size=1:100, shape=1, nsample=1e6, seed=0) {

  # bias for each number of successes
  bias = NULL

  # loop through numbers of successes
  for(i in size) {

    # random number seed
    set.seed(seed)

    # find the maximum bias over the success probabilities
    bias = c(bias, optimise(function(p) {

      # observed data
      PML.obs = 0
      k.obs = rnbinom(nsample, i, p)
      PML.obs = lchoose(k.obs+i-1, i-1) + lbeta(2*i+shape, 2*k.obs+shape)

      # replicate data
      k.rep = rnbinom(nsample, i, p)
      PML.rep = lchoose(k.rep+i-1, i-1) + lbeta(2*i+shape, k.rep+k.obs+shape)

      # the bias
      mean(PML.obs - PML.rep)

    }, c(0, 1), maximum=TRUE)$objective)
  }

  bias
}
