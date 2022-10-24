#' Compute EBF bias for F tests
#'
#' Computes the bias in the log posterior marginal likelihood for F tests.
#'
#' Pre-computed results are in the \code{ebf::f.bias} object,
#' for degrees of freedom from 1 to 100 in either numerator or denominator.
#' The bias is computed by numerically integrating over a bivariate F distriubtion.
#' This requires the \code{cubature} library.
#' If \code{df1} or \code{df2} is 1, the bias from the equivalent t-test is used.
#'
#' @param df1 Vector of numerator degrees of freedom for which to compute bias.
#'
#' @param df2 Vector of denominator degrees of freedom for which to compute bias.
#
#' @import stats
#'
#' @return A matrix with the biases corresponding to the elements of
#' \code{df1} and \code{df2}.  Numerator degrees of freedom are in rows, denominator degrees of
#' freedom in columns.
#'
#' @template references
#'
#' @export

compute.f.bias <- function(df1=1:100, df2=1:100) {

  # bias for each df
  bias = NULL

  # loop through sample sizes
  for(i in df1) {
    bias.tmp = NULL
    if (i>1)
      for(j in df2) {
        if (j>1)
          bias.tmp = c(bias.tmp, cubature::adaptIntegrate(function(x)
            df(x[1], i, j)* df(x[2], i, j) *
              (log(i/(i-1) * beta(i, j) / beta(i/2, j/2)^2 / x[1]) -
                 log(integrate(function(y)
                   df(x[1]*y, i, j) * df(x[2]*y, i, j),
                   0, Inf, stop.on.error=F)$value))
            , c(0,0), c(Inf, Inf))$integral)

        else bias.tmp = c(bias.tmp, compute.t.bias(i))
      }
    else bias.tmp = compute.t.bias(df2)
    bias = rbind(bias, bias.tmp)
  }
  rownames(bias) = NULL
  bias
}
