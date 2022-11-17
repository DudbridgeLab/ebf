#' Empirical Bayes factors for F tests
#'
#' Calculates empirical Bayes factors (EBFs) for F tests.
#'
#' The EBF includes bias adjustments to the log posterior marginal likelihoods.
#' Pre-computed adjustments are used for \code{df1} and \code{df2} from 1 to 100.
#' If either are 1, an equivalent two-sided t test is used.
#' For larger values, the asymptotic adjustment of 0.5 is used.
#' This may however be inaccurate when there is a large difference between
#' \code{df1} and \code{df2}.
#' In that case, \code{\link{compute.f.bias}} can be used to calculate an
#' accurate adjustment.  This can then be supplied here in the \code{bias} parameter.
#'
#' @template allParams
#'
#' @param df1 Vector of numerator degrees of freedom.
#'
#' @param df2 Vector of denominator degrees of freedom.
#'
#' @param bias Vector of biases for each test.
#'
#' @return A data frame containing the following components:
#'   \itemize{
#'   \item{\code{index} Indices of elements in \code{x} for which EBFs are calculated.}
#'     \item{\code{ebf} Empirical Bayes factors for the elements of \code{x}.}
#'     \item{\code{ebf.units} Expressed as units of evidence.}
#'     \item{\code{p} P-values.}
#'     \item{\code{p.log10} Expressed in -log10 units.}
#'   }
#'
#' @import stats
#'
#' @template references
#'
#' @export


ebf.f <- function(x,
                  df1, # no default value
                  df2,
                  index=NULL,
                  bias=NULL) {

  if (is.null(index)) index=1:length(x)
  df1 = df1 * rep(1, length(x))
  df2 = df2 * rep(1, length(x))

  ebf = rep(0, length(x))

  # posterior marginal likelihood
  for(i in 1:length(x)) {

    if (df1[i]>1 | df2[i]>1) {
      ebf[i] = exp(
        lbeta(df1[i], df2[i]) - 2*lbeta(df1[i]/2, df2[i]/2) +
          pf(x[i], 2*df1[i], 2*df2[i], log=T) - pf(x[i], df1[i], df2[i], log=T) -
          log(x[i]) -
          df(x[i], df1[i], df2[i], log=T)
      )
    }

  # bias
  if (is.null(bias))
    for(i in 1:length(x)) {
      if (df1[i] <= nrow(f.bias) & df2[i] <= ncol(f.bias))
        bias[i] = f.bias[df1[i], df2[i]]
      else
        bias[i] = compute.f.bias(df1[i], df2[i])
    }

  }

  # EBF
  ebf = ebf/exp(bias)
  ebf.units = log(ebf) / log((sqrt(3)+1)/(sqrt(3)-1))

  # P-value
  p = pf(x, df1, df2, lower=F)
  p.log10 = -log(p)/log(10)

  result = data.frame(index =  index,
                      ebf = ebf[index],
                      ebf.units = ebf.units[index],
                      p = p[index],
                      p.log10 = p.log10[index])

  result
}
