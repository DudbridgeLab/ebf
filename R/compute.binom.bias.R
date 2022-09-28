#' Compute EBF bias for binomial tests
#'
#' Computes the bias in the log posterior marginal likelihood for binomial tests.
#'
#' Pre-computed results are in the \code{ebf::binom.bias} object,
#' for sample sizes 1 to 100 with uniform prior.
#'
#' @param n Vector of sample sizes for which to estimate bias.
#' @param shape Parameter of symmetrical beta prior distribution.
#'
#' @import stats
#'
#' @return A data frame with the biases corresponding to the elements of
#' \code{n}.
#'
#' @author Frank Dudbridge
#'
#' @references
#' Dudbridge F (submitted) Units of evidence and expected Bayes factors for
#' objective reporting of statistical evidence.
#'
#' @export

compute.binom.bias <- function(n=1:100, shape=1) {

  # bias for each sample size
  bias = rep(0, length(n))

  # loop through sample sizes
  for(i in 1:length(n)) {

    # find the maximum bias over the success probabilities
    bias[i] = optimise(function(p) {
      # observed data
      PML.obs = 0
      PML.var = 0
      for(k.obs in 0:n[i]) {

        # replicate data
        PML.rep = 0
        for(k.rep in 0:n[i]) {
          # posterior marginal likelihood
          PML.rep = PML.rep + dbinom(k.rep, n[i], p) *
            (lchoose(n[i], k.rep) + lbeta(k.obs+k.rep+shape, 2*n[i]-k.obs-k.rep+shape))
        }

        # the bias for this observed data
        PML.obs = PML.obs + dbinom(k.obs, n[i], p) *
          (lchoose(n[i], k.obs) + lbeta(2*k.obs+shape, 2*(n[i]-k.obs)+shape) -
             PML.rep)
      }

      # expectation over observed data
      PML.obs
     }, c(0, 1), maximum=TRUE)$objective
  }

  data.frame(bias=bias, row.names=n)
}
