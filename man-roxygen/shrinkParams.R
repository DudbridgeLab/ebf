#' @param h0 If a scalar, the value of a point null hypothesis.
#' If a vector with two elements, the lower and upper
#' bounds of the null hypothesis.
#'
#' @param shrink If \code{TRUE}, shrinks EBFs to a posterior
#' distribution estimated from all the elements of \code{x}.
#'
#' @param points Vector of indices in \code{x} to use in estimating shrinkage EBFs.
#'
#' @param pi0 Proportion of true null hypotheses in shrinkage EBFs.
#'
#' @details
#' For shrinkage EBFs, each test is evaluated using the posterior distribution
#' obtained from every other test.
#' The calculation may be time-consuming when the number of tests is large.
#' To reduce the computation a subset of tests can be used instead,
#' with indices specified in \code{points}.  Each individual test is appended
#' to this list as it is evaluated.
#'
#' Further reduction can be achieved by calculating EBFs only for the elements
#' of \code{x} given by \code{index}.  The posterior distribution is still calculated using all
#' the points in \code{x}, or those specified in \code{points}.
#'
#' For test \code{i}, the contributions from tests \code{j!=i} can be weighted
#' by \code{pi0} to model the probability that the hypothesis is true for those
#' other tests.  The default \code{pi0=0} assumes that the hypothesis is true
#' for all tests.  If \code{pi0=1}, the single-test EBF is recovered.
#'
#' @return A data frame containing the following components:
#'   \itemize{
#'   \item{\code{index} Indices of elements in \code{x} for which EBFs are calculated.}
#'     \item{\code{ebf} Empirical Bayes factors for the elements of \code{x}.}
#'     \item{\code{ebf.units} Expressed as units of evidence, that is
#'      on the log scale with base 3.73.  One unit updates
#'      weaker belief to stronger belief, or weaker disbelief to weaker belief,
#'      according to a rational definition of weaker and stronger belief.
#'      See references for further details.}
#'     \item{\code{p} P-values if \code{h0} is a point null hypothesis.}
#'     \item{\code{p.log10} Expressed in -log10 units.}
#'     \item{\code{ebf.shrink} Shrunken EBFs if specified by \code{shrink==TRUE}.}
#'     \item{\code{ebf.shrink.units} Expressed as units of evidence.}
#'   }
