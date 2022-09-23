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
#' If \code{shrink==TRUE} the calculation is of order \code{length(x)^2}
#' and may be very time-consuming.  To reduce the computation a subset of
#' \code{x} may be used to calculate the shrinkage EBFs.
#' By default a maximum of 1000 points is used, taken from the quantiles of the
#' simple EBFs.  This number can be varied with the \code{npoints} option.
#'
#' It can be restricted to a subset of \code{x} by using \code{index}.
#' The shrunken EBFs are only calculated for the specified elements of \code{x},
#' however the full distribution of \code{x} is used in calculating those EBFs.
#'
#' Units of evidence are EBFs on the log scale with base 3.73.  One unit updates
#' weaker belief to stronger belief, or weaker disbelief to weaker belief,
#' according to an objective definition of weaker and stronger belief.
#' See Dudbridge (submitted) for further details.
#'
#' @template sharedParams
#'
#' @param test Type of EBF to compute.  Can be "normal" (default), "t" or "chisq".
#' @param df Vector of degrees of freedom for t or chi-squared tests.
#' Defaults to \code{Inf} for t test and 1 for chi-squared tests.
#' For the t test, values of \code{df} over 100 are set to \code{Inf}, equivalent to the Normal test.
#'
#' @import stats
#'
#' @template references
#'
#' @export

ebf <- function(x,
                se=NULL,
                h0=c(0,0),
                h1=NULL,
                test=c("normal","t","chisq"),
                df=NULL,
                shrink=FALSE,
                npoints=1000,
                index=NULL) {

  if (is.null(se)) se=1
  if (is.null(index)) index=1:length(x)
  test = tolower(test[1])

  if (startsWith("normal",test)) ebf = ebf.norm(x, se, h0, h1, shrink, npoints)

  if (startsWith("t",test)) ebf = ebf.t(x, se, h0, h1, df, shrink, npoints)

  if (startsWith("chisq",test)) ebf = ebf.chisq(x,df)

  ebf

}
