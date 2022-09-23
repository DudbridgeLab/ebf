#' Empirical Bayes factors for t tests
#'
#' Calculates empirical Bayes factors (EBFs) for univariate t tests.
#'
#' @template sharedParams
#'
#' @param df Vector of degrees of freedom.  Defaults to \code{Inf}.
#'
#' @import stats
#'
#' @template references
#'
#' @export

ebf.t <- function(x,
                  se=NULL,
                  h0=c(0,0),
                  h1=NULL,
                  df=NULL,
                  shrink=FALSE,
                  index=NULL,
                  npoints=1000) {

  if (is.null(df)) df = Inf

    # null hypothesis
  ### point hypothesis
  if (h0[1] == h0[2]) ebf.h0 = dt((x - h0[1]) /se, df) / se

    ### interval hypothesis
  if (h0[1] != h0[2]) ebf.h0 = ebf.t.simple(x, se, min(h0), max(h0), df)

  # alternative hypothesis
  if (!is.null(h1)) {
    ### point hypothesis
    if (h1[1] == h1[2]) ebf.h1 = dt((x - h1[1]) / se, df) / se
    ### interval hypothesis
    if (h1[1] != h1[2]) ebf.h1 = ebf.t.simple(x, se, min(h1), max(h1), df)
  }
  ### complement interval
  if (is.null(h1)) ebf.h1 = ebf.t.simple(x, se, min(h0), max(h0), df, TRUE)

    # EBFs
  ebf = ebf.h1 / ebf.h0
  ebf.units = log(ebf) / log((sqrt(3)+1)/(sqrt(3)-1))

  # P-values
  p = NULL
  if (h0[1]==0 & h0[2]==0) {
    if (!is.null(h1)) {
      ### two-sided test
      if (h1[1]==-Inf & h1[2]==Inf) p = pt(-abs(x) / se, df) *2
      ### one-sided positive test
      if (h1[1]==0 & h1[2]==Inf) p = pt(x / se, df, lower=F)
      ### one-sided negative test
      if (h1[1]==-Inf & h1[2]==0) p = pt(x / se, df)
    } else {
      ### two-sided test
      p = pt(-abs(x) / se, df) *2
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
      ebf.h0.shrink = ebf.t.shrink(x, se, min(h0), max(h0), points, df, index)

    # alternative hypothesis
    if (!is.null(h1)) {
      ### point hypothesis
      if (h1[1] == h1[2]) ebf.h1.shrink = ebf.h1[index]
      ### interval hypothesis
      if (h1[1] != h1[2])
        ebf.h1.shrink = ebf.t.shrink(x, se, min(h1), max(h1), df, points, index)
    }
    ### complement interval
    if (is.null(h1))
      ebf.h1.shrink = ebf.t.shrink(x, se, min(h0), max(h0), df, points, index, TRUE)

        ebf.shrink = ebf.h1.shrink / ebf.h0.shrink
    ebf.shrink.units = log(ebf.shrink) / log((sqrt(3)+1)/(sqrt(3)-1))
  }

  result = list(index = index,
                ebf = ebf,
                ebf.units = ebf.units,
                ebf.shrink = ebf.shrink,
                ebf.shrink.units = ebf.shrink.units,
                p = p,
                p.log10 = p.log10)
  class(result) = ("ebf")
  result
}
