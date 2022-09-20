#' Empirical Bayes factors
#'
#' Calculates empirical Bayes factors (EBFs) for a set of statistics.
#' Normal, t and chi-squared tests are available.
#' Two-sided, one-sided and interval hypotheses can be tested.
#'
#' This is the main function in the package.
#' By default, the statistics are normally distributed with known variances.
#' If \code{se} is unspecified then the variances are assumed to be 1.
#'
#' The default test has mean=0 as the null hypothesis, with two-sided alternative.
#' For one-sided alternatives, use \code{h1=c(0,-Inf)} or \code{h1=c(-Inf,0)}.
#' To test positive values against negative values, use \code{h0=c(-Inf,0)}.
#' By default \code{h1} is the complement of \code{h0}, in this case \code{c(0,Inf)}.
#' To test negative values against positive values, use \code{h0=c(0,Inf)}.
#'
#'
#' For chi-squared tests, \code{h0}, \code{h1} and \code{shrink} are ignored,
#'
#' If \code{shrink==TRUE} the calculation may be very
#' time-consuming when \code{length(x)} is large.
#' It can be restricted to a subset of \code{x} by using \code{index}.
#' The shrunken EBFs are only calculated for the specified elements of \code{x},
#' however the full set of statistics in \code{x} is used in calculating those EBFs.
#'
#' Units of evidence are EBFs on the log scale with base 3.73.  One unit updates
#' weaker belief to stronger belief, or weaker disbelief to weaker belief,
#' according to an objective definition of weaker and stronger belief.
#' See Dudbridge (submitted) for further details.
#'

#' @param x Vector of statistics for which EBFs will be calculated.
#' @param se Vector of standard errors associated with \code{x}.
#' @param h0 A vector with two elements, giving the lower and upper
#' bounds of the null hypothesis.
#' @param h1 A vector with two elements, giving the lower and upper bounds of
#' the alternative hypothesis to \code{h0}.  If unspecified, the alternative
#' hypothesis is the complement of \code{h0}.
#' @param test Type of EBF to compute.  Can be "Normal" (default), "t" or "chisq".
#' @param df Vector of degrees of freedom for t or chi-squared tests.
#' Defaults to \code{Inf} for t test and 1 for chi-squared tests.
#' For the t test, values of \code{df} over 30 are set to \code{Inf}, equivalent to the Normal test.
#' @param shrink If \code{TRUE}, uses information from all elements of \code{x}
#' when calculating each individual EBF.
#' @param index Vector of indices selecting a subset of \code{x}.
#'
#' @import stats
#'
#' @return An object of class "ebf" containing the following components:
#'   \itemize{
#'     \item{\code{ebf} Empirical Bayes factors for the elements of \code{x}.}
#'     \item{\code{ebf.units} \code{EBF} expressed as units of evidence.}
#'     \item{\code{ebf.shrink} Shrunken EBFs if specified by \code{shrink==TRUE}.}
#'     \item{\code{ebf.shrink.units} \code{ebf.shrink} expressed as units of evidence.}
#'     \item{\code{p} P-values.  This is empty if \code{h0} is an interval.}
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

ebf <- function(x,
                se=NULL,
                h0=c(0,0),
                h1=NULL,
                test=c("normal","t","chisq"),
                df=NULL,
                shrink=FALSE,
                index=NULL) {

  if (is.null(se)) se=1
  if (is.null(index)) index=1:length(x)
  test = tolower(test[1])

  if (startsWith("normal",test)) ebf = ebf.norm(x,se,h0,h1,shrink,index)

  if (startsWith("t",test)) ebf = ebf.t(x,se,h0,h1,df,shrink,index)

  if (startsWith("chisq",test)) ebf = ebf.chisq(x,df,index)

  ebf

}
