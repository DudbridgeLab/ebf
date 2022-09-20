#' Empirical Bayes factors for t tests
#'
#' Calculates empirical Bayes factors (EBFs) for univariate t tests.
#'
#' @param x Vector of statistics for which EBFs will be calculated.
#' @param se Vector of standard errors associated with \code{x}.
#' @param h0 A vector with two elements, giving the lower and upper
#' bounds of the null hypothesis.
#' @param h1 A vector with two elements, giving the lower and upper bounds of
#' the alternative hypothesis to \code{h0}.  If unspecified, the alternative
#' hypothesis is the complement of \code{h0}.
#' @param df Vector of degrees of freedom.  Defaults to \code{Inf}.
#' @param shrink If \code{TRUE}, uses information from all elements of \code{x}
#' when calculating each individual EBF.
#' @param index Vector of indices selecting a subset of \code{x}.

#' @import stats
#'
#' @return An object of class "ebf" containing the following components:
#'   \itemize{
#'     \item{\code{ebf} Empirical Bayes factors for the elements of \code{x}.}
#'     \item{\code{ebf.units} \code{ebf} expressed as units of evidence.}
#'     \item{\code{ebf.shrink} Shrunken EBFs if specified by \code{shrink==TRUE}.}
#'     \item{\code{ebf.shrink.units} \code{ebf.shrink} expressed as units of evidence.}
#'     \item{\code{p} P-values.}
#'     \item{\code{p.log10} -log10 of \code{p}.}
#'   }
#'
#' @author Frank Dudbridge
#'
#' @references
#' Dudbridge F (submitted) Units of evidence and expected Bayes factors for
#' objective reporting of statistical evidence.
#'
#' @export

ebf.t <- function(x,
                  se=NULL,
                  h0=c(0,0),
                  h1=NULL,
                  df=NULL,
                  shrink=FALSE,
                  index=NULL) {
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
  ebf.h0 = ebf.h0[index]
  ebf.h1 = ebf.h1[index]
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
  p = p[index]
  p.log10 = NULL
  if (!is.null(p)) p.log10 = -log(p)/log(10)

  # shrinkage
  ebf.shrink = NULL
  ebf.shrink.units = NULL
  if (shrink == TRUE) {

    se = rep(0,length(x)) + se # make into full length vector

    # null hypothesis
    ### point hypothesis
    if (h0[1] == h0[2]) ebf.h0.shrink = ebf.h0
    ### interval hypothesis
    if (h0[1] != h0[2])
      ebf.h0.shrink = ebf.t.shrink(x, se, min(h0), max(h0), df, index)

    # alternative hypothesis
    if (!is.null(h1)) {
      ### point hypothesis
      if (h1[1] == h1[2]) ebf.h1.shrink = ebf.h1
      ### interval hypothesis
      if (h1[1] != h1[2])
        ebf.h1.shrink = ebf.t.shrink(x, se, min(h1), max(h1), df, index)
    }
    ### complement interval
    if (is.null(h1))
      ebf.h1.shrink = ebf.t.shrink(x, se, min(h0), max(h0), df, index, TRUE)

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
