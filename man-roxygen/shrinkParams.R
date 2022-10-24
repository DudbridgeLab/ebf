#' @param h0 If a scalar, the value of a point null hypothesis.
#' If a vector with two elements, the lower and upper
#' bounds of the null hypothesis.
#'
#' @param shrink If \code{TRUE}, shrinks EBFs to a posterior
#' distribution estimated from all the elements of \code{x}.
#'
#' @param npoints Number of elements of \code{x} to use in estimating posterior distribution.
#'
#' @param nsupport Number of mixture components in posterior distribution.
#'
#' @param tol Convergence tolerance in EM algorithm for fitting posterior distribution.
#'
#' @param nboot Number of samples in bootstrap estimation of posterior distribution.
#'
#' @param seed Random number seed.
#'
#' @details
#' For shrinkage EBFs, two approaches are available for estimating the posterior
#' distribution.  The default, generally faster option fits a mixture of conjugate
#' distributions to \code{x}. All components of the mixture are given the same fixed precision,
#' which in general terms is the precision of the posterior when a single component
#' is fitted.  (For the t-tests, which do not have a conjugate distribution, a mixture
#' of t distributions is fitted, reflecting a flat prior.)
#' The number of mixture components can be varied by \code{nsupport}.
#'
#' If \code{nboot > 0}, a discrete distribution is fitted, corresponding to the
#' maximum likelihood non-parametric estimate of the distribution of \code{x}.
#' This is also equivalent to the mixture of conjugate distributions where the
#' precision is fixed to be infinite.  To account for variation in this estimate,
#' EBFs are averaged over parametric bootstrap resamples of \code{x}.
#' The number of support points in the fitted distribution can be varied by
#' \code{nsupport}.
#'
#' The calculation may be time-consuming.  To reduce the computation a random subset of
#' \code{x} may be used to estimate the posterior distribution.
#' By default a maximum of 1000 points is used.
#' This number can be varied with the \code{npoints} parameter.
#'
#' Further reduction can be achieved by calculating EBFs only for the elements
#' of \code{x} given by \code{index}.  Those EBFs are calculated using all
#' the points in \code{x}, unless varied by \code{npoints}.
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
