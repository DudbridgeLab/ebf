#' Empirical Bayes factors for Poisson tests
#'
#' Calculates empirical Bayes factors (EBFs) for Poisson tests of rates.
#'
#' The EBF includes bias adjustments to the log posterior marginal likelihoods.
#' The adjustment depends on the unknown rate parameter.
#' The maximum value is 0.5346, corresponding to a rate of 7.301.
#' This is used to give a conservative EBF.
#'
#' If a normal approximation is acceptable, use \code{\link{ebf.norm}}.
#'
#' @template allParams
#' @template shrinkParams
#'
#' @param interval Vector containing the interval lengths in each test.
#' The numbers of events are contained in \code{x}.
#'
#' @param h1 If a scalar, the value of a point alternative hypothesis.
#' If a vector with two elements, the lower and upper bounds of
#' the alternative hypothesis.  If \code{NULL} (default), the alternative
#' hypothesis is the complement of \code{h0}.
#'
#' The default test has \code{h0=1}, with two-sided alternative.
#' For one-sided alternatives, use \code{h1=c(1,Inf)} or \code{h1=c(0,1)}.
#' To test higher values against lower values, use \code{h0=c(0,1)}.
#' In this case \code{h1} defaults to \code{c(1,Inf)}.
#' To test lower values against higher values, use \code{h0=c(1,Inf)}.
#'
#'
#' @import stats
#'
#' @template references
#'
#' @export

ebf.poisson <- function(x,
                  interval=1,
                  index=NULL,
                  h0=1,
                  h1=NULL,
                  shrink=FALSE,
                  npoints=1000
                  ) {

  if (length(h0) == 1) h0 = c(h0, h0)
  if (length(h1) == 1) h1 = c(h1, h1)
  if (is.null(index)) index=1:length(x)

  # expand to a vector
  interval = rep(0, length(x)) + interval

  # null hypothesis
  ### point hypothesis
  if (h0[1] == h0[2]) ebf.h0 = dpois(x, interval * h0[1])

  ### interval hypothesis
  if (h0[1] != h0[2]) ebf.h0 = ebf.poisson.simple(x, interval, min(h0), max(h0))

  # alternative hypothesis
  if (!is.null(h1)) {
    ### point hypothesis
    if (h1[1] == h1[2]) ebf.h1 = dpois(x, interval * h1[1])
    ### interval hypothesis
    if (h1[1] != h1[2]) ebf.h1 = ebf.poisson.simple(x, interval, min(h1), max(h1))
  }
  ### complement interval
  if (is.null(h1)) ebf.h1 = ebf.poisson.simple(x, interval, min(h0), max(h0), TRUE)

    # EBFs
  ebf = ebf.h1 / ebf.h0
  ebf.units = log(ebf) / log((sqrt(3)+1)/(sqrt(3)-1))

  # P-values
  p = rep(NA,length(x))
  if (h0[1]==h0[2]) {
    if (!is.null(h1)) {
      ### two-sided test
      if (h1[1]==0 & h1[2]==Inf) {
        for(i in 1:length(x))
          p[i] = poisson.test(x[i], interval[i], h0[1])$p.value
      }
      ### one-sided higher test
      if (h1[1]==h0[1] & h1[2]==1) p = ppois(x-1, interval * h0[1], lower=F)
      ### one-sided negative test
      if (h1[1]==0 & h1[2]==h0[2]) p = ppois(x, interval * h0[1])
    } else {
      ### two-sided test
      for(i in 1:length(x))
        p[i] = poisson.test(x[i], interval[i], h0[1])$p.value
    }
  }
  p.log10 = NULL
  if (!is.null(p)) p.log10 = -log(p)/log(10)

  # shrinkage
  ebf.shrink = NULL
  ebf.shrink.units = NULL
  if (shrink == TRUE) {

    # select points to use in estimating shrinkage EBFs
    if (length(x) < npoints) points = 1:length(x)
    else points = order(ebf)[seq(1/npoints, 1, 1/npoints)*length(x)]

    # null hypothesis
    ### point hypothesis
    if (h0[1] == h0[2]) ebf.h0.shrink = ebf.h0[index]
    ### interval hypothesis
    if (h0[1] != h0[2])
      ebf.h0.shrink = ebf.poisson.shrink(x, interval, index, min(h0), max(h0), points)

    # alternative hypothesis
    if (!is.null(h1)) {
      ### point hypothesis
      if (h1[1] == h1[2]) ebf.h1.shrink = ebf.h1[index]
      ### interval hypothesis
      if (h1[1] != h1[2])
        ebf.h1.shrink = ebf.poisson.shrink(x, interval, index, min(h1), max(h1), points)
    }
    ### complement interval
    if (is.null(h1))
      ebf.h1.shrink = ebf.poisson.shrink(x, interval, index, min(h0), max(h0), points, TRUE)

    ebf.shrink = ebf.h1.shrink / ebf.h0.shrink
    ebf.shrink.units = log(ebf.shrink) / log((sqrt(3)+1)/(sqrt(3)-1))
  }

  result = data.frame(index =  index,
                      ebf = ebf[index],
                      ebf.units = ebf.units[index],
                      p = p[index],
                      p.log10 = p.log10[index])

  if (shrink == TRUE) result = data.frame(result,
                                          ebf.shrink,
                                          ebf.shrink.units)

  result
}
