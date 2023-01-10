#' Compute EBF bias for F tests
#'
#' Computes the bias in the log posterior marginal likelihood for F tests.
#'
#' Pre-computed results are in the \code{ebf::f.bias} object,
#' for degrees of freedom from 1 to 100 in either numerator or denominator.
#' The bias is computed by numerically integrating over a bivariate F distriubtion.
#' This requires the \code{cubature} library.
#'
#' @param df1 Vector of numerator degrees of freedom for which to compute bias.
#'
#' @param df2 Vector of denominator degrees of freedom for which to compute bias.
#
#' @import stats
#'
#' @return A matrix with the biases corresponding to the elements of
#' \code{df1} and \code{df2}.
#'
#' @template references
#'
#' @export

compute.f.bias <- function(df1=1:100, df2=1:100) {

  # bias for each df
  bias = matrix(nrow=length(df1), ncol=length(df2))

  # loop through sample sizes
  for(i in 1:length(df1)) {
    for(j in 1:length(df2)) {
      k1 = df1[i]
      k2 = df2[j]
      #if (k1+k2>2)
        bias[i,j] = lbeta(k1, k2) - 2*lbeta(k1/2, k2/2) +
        cubature::cubintegrate(function(x)
          -log(x) * df(x, k1, k2), 0, Inf)$integral -
        cubature::cubintegrate(function(x)
          df(x[1], k1, k2) * df(x[2], k1, k2) *
            log(
              cubature::hcubature(function(y)
              #integrate(function(y)
                x[1] * df(x[1]*y, k1, k2) * y * df(x[2]*y, k1, k2),
                0, Inf)$integral),
                #0, Inf, stop.on.error=TRUE)$value),
          c(0,0), c(Inf, Inf), absTol =1e-4)$integral
    }
  }

  bias
}
