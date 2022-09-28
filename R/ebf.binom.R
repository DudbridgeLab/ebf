#' Empirical Bayes factors for binomial tests
#'
#' Calculates empirical Bayes factors (EBFs) for binomial tests of proportions.
#'
#' This provides an exact but conservative EBF.
#' If a normal approximation to the binomial is acceptable,
#' \code{\link{ebf.norm}} is faster.
#'
#' @template allParams
#' @template shrinkParams
#'
#' @param size Vector containing the numbers of trials in each test.
#' The numbers of successes are contained in \code{x}.
#'
#' @param h1 If a scalar, the value of a point alternative hypothesis.
#' If a vector with two elements, the lower and upper bounds of
#' the alternative hypothesis.  If \code{NULL} (default), the alternative
#' hypothesis is the complement of \code{h0}.
#'
#' The default test has \code{h0=0.5}, with two-sided alternative.
#' For one-sided alternatives, use \code{h1=c(0.5,1)} or \code{h1=c(0,0.5)}.
#' To test higher values against lower values, use \code{h0=c(0,0.5)}.
#' In this case \code{h1} defaults to \code{c(0.5,1)}.
#' To test lower values against higher values, use \code{h0=c(0.5,1)}.
#'
#' @param shape Parameter of the prior distribution Beta(\code{shape},\code{shape}).
#'
#' @import stats
#'
#' @template references
#'
#' @export

ebf.binom <- function(x,
                  size,
                  index=NULL,
                  h0=0.5,
                  h1=NULL,
                  shrink=FALSE,
                  npoints=1000,
                  shape=1
                  ) {

  if (length(h0) == 1) h0 = c(h0, h0)
  if (length(h1) == 1) h1 = c(h1, h1)
  if (is.null(index)) index=1:length(x)

  # expand to a vector
  size = rep(0, length(x)) + size

  # null hypothesis
  ### point hypothesis
  if (h0[1] == h0[2]) ebf.h0 = dbinom(x, size, h0[1])

  ### interval hypothesis
  if (h0[1] != h0[2]) ebf.h0 = ebf.binom.simple(x, size, min(h0), max(h0), shape)

  # alternative hypothesis
  if (!is.null(h1)) {
    ### point hypothesis
    if (h1[1] == h1[2]) ebf.h1 = dbinom(x, size, h1[1])
    ### interval hypothesis
    if (h1[1] != h1[2]) ebf.h1 = ebf.binom.simple(x, size, min(h1), max(h1), shape)
  }
  ### complement interval
  if (is.null(h1)) ebf.h1 = ebf.binom.simple(x, size, min(h0), max(h0), shape, TRUE)

    # EBFs
  ebf = ebf.h1 / ebf.h0
  ebf.units = log(ebf) / log((sqrt(3)+1)/(sqrt(3)-1))

  # P-values
  p = NULL
  if (h0[1]==h0[2]) {
    if (!is.null(h1)) {
      ### two-sided test
      if (h1[1]==0 & h1[2]==1) p = pbinom(min(x, size-x), size, h0[1]) *2
      ### one-sided higher test
      if (h1[1]==h0[1] & h1[2]==1) p = pbinom(size-x, size, h0[1])
      ### one-sided negative test
      if (h1[1]==0 & h1[2]==h0[1]) p = pbinom(x, size, h0[1])
    } else {
      ### two-sided test
      p = pbinom(min(x, size-x), size, h0[1]) *2
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
      ebf.h0.shrink = ebf.binom.shrink(x, size, index, min(h0), max(h0), points, shape)

    # alternative hypothesis
    if (!is.null(h1)) {
      ### point hypothesis
      if (h1[1] == h1[2]) ebf.h1.shrink = ebf.h1[index]
      ### interval hypothesis
      if (h1[1] != h1[2])
        ebf.h1.shrink = ebf.binom.shrink(x, size, index, min(h1), max(h1), points, shape)
    }
    ### complement interval
    if (is.null(h1))
      ebf.h1.shrink = ebf.binom.shrink(x, size, index, min(h0), max(h0), points, shape, TRUE)

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
