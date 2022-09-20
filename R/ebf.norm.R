#' Empirical Bayes factors for normal tests
#'
#' Calculates empirical Bayes factors (EBFs) for univariate normal tests.
#'
#' @param x Vector of statistics for which EBFs will be calculated.
#' @param se Vector of standard errors associated with \code{x}.
#' @param h0 A vector with two elements, giving the lower and upper
#' bounds of the null hypothesis.
#' @param h1 A vector with two elements, giving the lower and upper bounds of
#' the alternative hypothesis to \code{h0}.  If unspecified, the alternative
#' hypothesis is the complement of \code{h0}.
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

ebf.norm <- function(x,
                     se=NULL,
                     h0=c(0,0),
                     h1=NULL,
                     shrink=FALSE,
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
  ebf.h0 = ebf.h0[index]
  ebf.h1 = ebf.h1[index]
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
    if (h0[1] == h0[2]) ebf.h0.shrink = dnorm(x, h0[1], se) [index]
    ### interval hypothesis
    if (h0[1] != h0[2])
      ebf.h0.shrink = ebf.norm.shrink(x, se, min(h0), max(h0), index)

    # alternative hypothesis
    if (!is.null(h1)) {
      ### point hypothesis
      if (h1[1] == h1[2]) ebf.h1.shrink = dnorm(x, h1[1], se) [index]
      ### interval hypothesis
      if (h1[1] != h1[2])
        ebf.h1.shrink = ebf.norm.shrink(x, se, min(h1), max(h1), index)
    }
    ### complement interval
    if (is.null(h1))
      ebf.h1.shrink = ebf.norm.shrink(x, se, min(h0), max(h0), index, TRUE)

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
