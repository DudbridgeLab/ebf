#' Empirical Bayes factors for Poisson tests
#'
#' Calculates empirical Bayes factors (EBFs) for one-sample Poisson tests of rates.
#'
#' @template allParams
#'
#' @param interval Vector containing the interval lengths in each test.
#' The numbers of events are contained in \code{x}.
#'
#' @param h1 If a scalar, the value of a point alternative hypothesis.
#' If a vector with two elements, the lower and upper bounds of
#' the alternative hypothesis.  If \code{NULL} (default), the alternative
#' hypothesis is the complement of \code{h0}.
#'
#' @param shape Shape parameter of the prior Gamma distribution.
#' @param rate Rate parameter of the prior Gamma distribution.
#'
#' @details
#' The default test has \code{h0=1}, with two-sided alternative.
#' For one-sided alternatives, use \code{h1=c(1,Inf)} or \code{h1=c(0,1)}.
#' To test higher values against lower values, use \code{h0=c(0,1)}.
#' In this case \code{h1} defaults to \code{c(1,Inf)}.
#' To test lower values against higher values, use \code{h0=c(1,Inf)}.
#'
#' @template shrinkParams
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
                        shape=1,
                        rate=0,
                        points=NULL,
                        pi0=0) {

  if (length(h0) == 1) h0 = c(h0, h0)
  if (length(h1) == 1) h1 = c(h1, h1)
  if (is.null(index)) index=1:length(x)

  # expand to a vector
  interval = rep(0, length(x)) + interval

  # null hypothesis
  ### point hypothesis
  if (h0[1] == h0[2]) ebf.h0 = dpois(x, interval * h0[1])

  ### interval hypothesis
  if (h0[1] != h0[2]) ebf.h0 = ebf.poisson.simple(x, interval,
                                                  min(h0), max(h0), shape, rate)

  # alternative hypothesis
  if (!is.null(h1)) {
    ### point hypothesis
    if (h1[1] == h1[2]) ebf.h1 = dpois(x, interval * h1[1])
    ### interval hypothesis
    if (h1[1] != h1[2]) ebf.h1 = ebf.poisson.simple(x, interval,
                                                    min(h1), max(h1), shape, rate)
  }
  ### complement interval
  if (is.null(h1)) ebf.h1 = ebf.poisson.simple(x, interval,
                                               min(h0), max(h0), shape, rate, TRUE)

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

    # data points for estimating non-parametric distribution
    if (is.null(points)) points = 1:length(x)

    # null hypothesis
    ### point hypothesis
    if (h0[1] == h0[2]) ebf.h0.shrink = ebf.h0
    ### interval hypothesis
    if (h0[1] != h0[2])
      ebf.h0.shrink = ebf.poisson.shrink(x, interval, index, min(h0), max(h0),
                                         shape, rate, points)

    # alternative hypothesis
    if (!is.null(h1)) {
      ### point hypothesis
      if (h1[1] == h1[2]) ebf.h1.shrink = ebf.h1
      ### interval hypothesis
      if (h1[1] != h1[2]) {
        if (h0[1] == h0[2])
          ebf.h1.shrink = ebf.poisson.shrink(x, interval, index, min(h1), max(h1),
                                             shape, rate, points, pi0)
        else
          ebf.h1.shrink = ebf.poisson.shrink(x, interval, index, min(h1), max(h1),
                                             shape, rate, points)
      }
    }
    ### complement interval
    if (is.null(h1)) {
      if (h0[1] == h0[2])
        ebf.h1.shrink = ebf.poisson.shrink(x, interval, index, min(h0), max(h0),
                                           shape, rate, points, pi0, TRUE)
      else
        ebf.h1.shrink = ebf.poisson.shrink(x, interval, index, min(h0), max(h0),
                                           shape, rate, points, complement=TRUE)
    }

    ebf.shrink = ebf.h1.shrink / ebf.h0.shrink
    ebf.shrink.units = log(ebf.shrink) / log((sqrt(3)+1)/(sqrt(3)-1))
  }

  result = data.frame(index =  index,
                      ebf = ebf[index],
                      ebf.units = ebf.units[index])

  if (sum(!is.na(p[index])) > 0) result = data.frame(result,
                                                     p = p[index],
                                                     p.log10 = p.log10[index])

  if (shrink == TRUE) result = data.frame(result,
                                          ebf.shrink = ebf.shrink[index],
                                          ebf.shrink.units = ebf.shrink.units[index])

  result
}
