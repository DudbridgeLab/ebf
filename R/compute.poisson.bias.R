#' Compute EBF bias for Poisson tests
#'
#' Computes the bias in the log posterior marginal likelihood for Poisson tests.
#'
#' @import stats
#'
#' @return The maximum bias over the rate parameter, which is 0.5003284.
#'
#' @template references
#'
#' @export

compute.poisson.bias <- function() {
  optimise(function(l) {
    bias=0
    for(i in 0:1000)
      for(j in 0:1000) {
        bias = bias + dpois(i, l) * dpois(j, l) *
          ((j-i)*log(2) + lfactorial(2*i) + lfactorial(j+1) - lfactorial(i+j) - lfactorial(i+1))
      }
    bias
  }, c(0,20), maximum=T)#$objective
}
