#' Empirical Bayes factors for binomial tests
#'
#' Calculates empirical Bayes factors (EBFs) for one-sample binomial tests of proportions.
#'
#' @template allParams
#' @template binomParams
#' @template shrinkParams
#'
#' @param size Vector containing the numbers of trials in each test.
#' The numbers of successes are contained in \code{x}.
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
                      shape=1,
                      points=NULL,
                      pi0=0) {

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
  p = rep(NA, length(x))
  if (h0[1]==h0[2]) {
    if (!is.null(h1)) {
      ### two-sided test
      if (h1[1]==0 & h1[2]==1) {
        for(i in 1:length(x))
          p[i] = as.numeric(binom.test(x[i], size[i], h0[1])$p.value)
      }
      ### one-sided higher test
      if (h1[1]==h0[1] & h1[2]==1) p = as.numeric(pbinom(x-1, size, h0[1], lower=F))
      ### one-sided lower test
      if (h1[1]==0 & h1[2]==h0[2]) p = as.numeric(pbinom(x, size, h0[1]))
    } else {
      ### two-sided test
      for(i in 1:length(x))
        p[i] = as.numeric(binom.test(x[i], size[i], h0[1])$p.value)
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
      ebf.h0.shrink = ebf.binom.shrink(x, size, index, min(h0), max(h0), shape, FALSE,
                                     points, nsupport, tol, correct)

    # alternative hypothesis
    if (!is.null(h1)) {
      ### point hypothesis
      if (h1[1] == h1[2]) ebf.h1.shrink = ebf.h1
      ### interval hypothesis
      if (h1[1] != h1[2]) {
        if (h0[1] == h0[2])
          ebf.h1.shrink = ebf.binom.shrink(x, size, index, min(h1), max(h1), shape,
                                         points, pi0)
        else
          ebf.h1.shrink = ebf.binom.shrink(x, size, index, min(h1), max(h1), shape,
                                         points)
      }
    }
    ### complement interval
    if (is.null(h1)) {
      if (h0[1] == h0[2])
        ebf.h1.shrink = ebf.binom.shrink(x, size, index, min(h0), max(h0), shape,
                                       points, pi0, TRUE)
      else
        ebf.h1.shrink = ebf.binom.shrink(x, size, index, min(h0), max(h0), shape,
                                       points, complement=TRUE)
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
