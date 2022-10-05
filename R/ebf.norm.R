#' Empirical Bayes factors for normal tests
#'
#' Calculates empirical Bayes factors (EBFs) for univariate normal tests.
#'
#' @template allParams
#' @template normParams
#' @template shrinkParams
#'
#' @import stats
#'
#' @template references
#'
#' @export

ebf.norm <- function(x,
                     se=1,
                     index=NULL,
                     h0=0,
                     h1=NULL,
                     shrink=FALSE,
                     npoints=1000) {

  if (length(h0) == 1) h0 = c(h0, h0)
  if (length(h1) == 1) h1 = c(h1, h1)
  if (is.null(index)) index=1:length(x)

  # expand into a vector
  se = rep(0,length(x)) + se

  # null hypothesis
  ### point hypothesis
  if (h0[1] == h0[2]) ebf.h0 = dnorm(x, h0[1], se)
  ### interval hypothesis
  if (h0[1] != h0[2]) ebf.h0 = ebf.norm.simple(x, se, min(h0), max(h0))

  # alternative hypothesis
  if (!is.null(h1)) {
    ### point hypothesis
    if (h1[1] == h1[2]) ebf.h1 = dnorm(x, h1[1], se)
    ### interval hypothesis
    if (h1[1] != h1[2]) ebf.h1 = ebf.norm.simple(x, se, min(h1), max(h1))
  }
  ### complement interval
  if (is.null(h1)) ebf.h1 = ebf.norm.simple(x, se, min(h0), max(h0), TRUE)

  # EBFs
  ebf = ebf.h1 / ebf.h0
  ebf.units = log(ebf) / log((sqrt(3)+1)/(sqrt(3)-1))

  # P-values
  p = rep(NA, length(x))
  if (h0[1]==h0[2]) {
    if (!is.null(h1)) {
      ### two-sided test
      if (h1[1]==-Inf & h1[2]==Inf) p = pnorm(-abs(x-h0[1]), 0, se) *2
      ### one-sided positive test
      if (h1[1]==h0[1] & h1[2]==Inf) p = pnorm(x-h0[1], 0, se, lower=F)
      ### one-sided negative test
      if (h1[1]==-Inf & h1[2]==h0[2]) p = pnorm(x-h0[1], 0, se)
    } else {
      ### two-sided test
      p = pnorm(-abs(x-h0[1]), 0 , se) *2
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
      ebf.h0.shrink = ebf.norm.shrink(x, se, index, min(h0), max(h0), points)

    # alternative hypothesis
    if (!is.null(h1)) {
      ### point hypothesis
      if (h1[1] == h1[2]) ebf.h1.shrink = ebf.h1[index]
      ### interval hypothesis
      if (h1[1] != h1[2])
        ebf.h1.shrink = ebf.norm.shrink(x, se, index, min(h1), max(h1), points)
    }
    ### complement interval
    if (is.null(h1))
      ebf.h1.shrink = ebf.norm.shrink(x, se, index, min(h0), max(h0), points, TRUE)

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
