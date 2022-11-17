#' Compute EBF bias for t tests
#'
#' Computes the bias in the log posterior marginal likelihood for t tests.
#'
#' Pre-computed results for 1 to 100 degrees of freedom are in the
#' \code{ebf::t.bias} object.
#' The bias is computed by numerically integrating over a bivariate t distriubtion.
#' This requires the \code{cubature} library.
#'
#' @param df Vector of degrees of freedom for which to compute bias.
#'
#' @import stats
#'
#' @return A vector with the biases corresponding to \code{df}.
#'
#' @template references
#'
#' @export

compute.t.bias <- function(df=1:100) {

  # posterior marginal likelihood, log scale
  PBF = 2*lgamma((df+1)/2) + lgamma(df+0.5) -
    (2*lgamma(df/2) + lgamma(df+1) + log(df*pi)/2)

  # expected marginal likelihood, over observed and replicate data, log scale
  EBF = sapply(df, function(k) {
    log.k.pi = log(k*pi)
    cubature::adaptIntegrate(function(x) { # integrate over obs and rep data
      a = ((1+x[1]^2 / k) * (1+x[2]^2 / k)) ^ (-k/2-0.5) # t density unnormalised
      if(a > 0)
        b = #integrate(function(z)
          cubature::cubintegrate(function(z)
          ((1+(x[1]-z)^2 / k) * (1+(x[2]-z)^2 / k))^(-k/2-0.5), # marginal likelihood
          #-Inf, Inf, stop.on.error=TRUE)$value
          -Inf, Inf)$integral
      else b = 0
      if (b > 0) a*(log(b) + 2*(lgamma(k/2+0.5) - lgamma(k/2)) - log.k.pi) # normalised
      else 0
    }, c(-Inf, -Inf), c(Inf, Inf))$integral

  })
  EBF = EBF * ( gamma((df+1)/2) / gamma(df/2) / sqrt(df*pi) )^2 # normalised

  bias = PBF - EBF

  bias

}
