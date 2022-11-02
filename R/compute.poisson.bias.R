#' Compute EBF bias for Poisson tests
#'
#' Computes the bias in the log posterior marginal likelihood for Poisson tests.
#'
#' @import stats
#'
#' @return The maximum bias over the rate parameter, which is 0.541.
#'
#' @template references
#'
#' @export

compute.poisson.bias <- function() {
  optimise(function(l) {
    bias=0
    for(i in 1:100)
      for(j in 0:100) {
        bias = bias + dpois(i, l)/(1-dpois(0, l)) * dpois(j, l) *
          ((j-i)*log(2) + lgamma(2*i) + lgamma(j+1) - lgamma(i+j) - lgamma(i+1))
      }
    bias
  }, c(0,20), maximum=T)$objective
}
