#' Compute EBF bias for binomial tests
#'
#' Computes the bias in the log posterior marginal likelihood for binomial tests.
#'
#' Pre-computed results are in the \code{ebf::binom.bias} object,
#' for sample sizes 1 to 100 with uniform prior.
#'
#' @param n Vector of sample sizes for which to compute bias.
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

compute.binom.bias <- function(n=1:100, shape=1) {

  # bias for each sample size
  bias = NULL

  # loop through sample sizes
  for(i in n) {

    # find the maximum bias over the success probabilities
    #bias = c(bias, optimise(function(p) {
    bias = c(bias, integrate(function(p) {
      # observed data
      PML.obs = 0
      PML.var = 0
      for(k.obs in 0:i) {

        # replicate data
        PML.rep = 0
        for(k.rep in 0:i) {
          # posterior marginal likelihood
          PML.rep = PML.rep + dbinom(k.rep, i, p) *
            #(lchoose(i, k.rep) + lbeta(k.obs+k.rep+shape, 2*i-k.obs-k.rep+shape))
            log(integrate(function(pp) dbinom(k.rep, i, pp) * dbeta(pp, k.obs+shape, i-k.obs+shape),
                      0.3, 0.7)$value)
        }

        # the bias for this observed data
        PML.obs = PML.obs + dbinom(k.obs, i, p) *
#          (lchoose(i, k.obs) + lbeta(2*k.obs+shape, 2*(i-k.obs)+shape) -
#             PML.rep)
        (log(integrate(function(pp) dbinom(k.obs, i, pp) * dbeta(pp, k.obs+shape, i-k.obs+shape),
                       0.3, 0.7)$value) - PML.rep)
        PML.var = PML.var + dbinom(k.obs, i, p) *
          (lchoose(i, k.obs) + lbeta(2*k.obs+shape, 2*(i-k.obs)+shape) -
             PML.rep)^2
      }

      # expectation over observed data
      PML.obs
      #PML.var - PML.obs^2
     #}, c(0, 1), maximum=TRUE)$objective)
    }, 0, 1)$value)
  }

  bias
}
