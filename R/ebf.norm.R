#' Empirical Bayes factors for normal tests
#'
#' Calculates empirical Bayes factors (EBFs) for univariate normal tests.
#'
#' @template sharedParams
#'
#' @import stats
#'
#' @template references
#'
#' @export

ebf.norm <- function(x,
                     se=NULL,
                     h0=c(0,0),
                     h1=NULL,
                     shrink=FALSE,
                     npoints=1000,
                     index=NULL) {

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
  p = NULL
  if (h0[1]==0 & h0[2]==0) {
    if (!is.null(h1)) {
      ### two-sided test
      if (h1[1]==-Inf & h1[2]==Inf) p = pnorm(-abs(x), 0, se) *2
      ### one-sided positive test
      if (h1[1]==0 & h1[2]==Inf) p = pnorm(x, 0, se, lower=F)
      ### one-sided negative test
      if (h1[1]==-Inf & h1[2]==0) p = pnorm(x, 0, se)
    } else {
      ### two-sided test
      p = pnorm(-abs(x), 0, se) *2
    }
  }
  p.log10 = NULL
  if (!is.null(p)) p.log10 = -log(p)/log(10)

  # shrinkage
  ebf.shrink = NULL
  ebf.shrink.units = NULL
  if (shrink == TRUE) {

    se = rep(0,length(x)) + se # make into full length vector
    if (is.null(index)) index = 1:length(x)

    # select points to use in estimating shrinkage EBFs
    if (length(ebf) < npoints) points = 1:length(ebf)
    else points = order(ebf)[seq(1/npoints, 1, 1/npoints)*length(ebf)]

    # null hypothesis
    ### point hypothesis
    if (h0[1] == h0[2]) ebf.h0.shrink = ebf.h0[index]
    ### interval hypothesis
    if (h0[1] != h0[2])
      ebf.h0.shrink = ebf.norm.shrink(x, se, min(h0), max(h0), points, index)

    # alternative hypothesis
    if (!is.null(h1)) {
      ### point hypothesis
      if (h1[1] == h1[2]) ebf.h1.shrink = ebf.h1[index]
      ### interval hypothesis
      if (h1[1] != h1[2])
        ebf.h1.shrink = ebf.norm.shrink(x, se, min(h1), max(h1), points, index)
    }
    ### complement interval
    if (is.null(h1))
      ebf.h1.shrink = ebf.norm.shrink(x, se, min(h0), max(h0), points, index, TRUE)

    ebf.shrink = ebf.h1.shrink / ebf.h0.shrink
    ebf.shrink.units = log(ebf.shrink) / log((sqrt(3)+1)/(sqrt(3)-1))
  }

  result = list(index =  index,
                ebf = ebf,
                ebf.units = ebf.units,
                ebf.shrink = ebf.shrink,
                ebf.shrink.units = ebf.shrink.units,
                p = p,
                p.log10 = p.log10)
  class(result) = ("ebf")
  result
}
