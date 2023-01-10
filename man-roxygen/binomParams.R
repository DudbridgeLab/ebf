#' @param h1 If a scalar, the value of a point alternative hypothesis.
#' If a vector with two elements, the lower and upper bounds of
#' the alternative hypothesis.  If \code{NULL} (default), the alternative
#' hypothesis is the complement of \code{h0}.
#'
#' @param shape Parameter of the symmetric prior Beta distribution.
#'
#' @details
#' The default test has \code{h0=0.5}, with two-sided alternative.
#' For one-sided alternatives, use \code{h1=c(0.5,1)} or \code{h1=c(0,0.5)}.
#' To test higher values against lower values, use \code{h0=c(0,0.5)},
#' in which case \code{h1} defaults to \code{c(0.5,1)}.
#' To test lower values against higher values, use \code{h0=c(0.5,1)}.
#'
#' The EBF includes bias adjustments to the log posterior marginal likelihoods.
#' When the hypothesis is \code{[0,1]}, pre-computed adjustments are used for \code{size} from 1 to 100.
#' For higher values of \code{size}, the asymptotic adjustment of 0.5 is used.
#'
#' If a normal approximation is acceptable, use \code{\link{ebf.norm}}.
#'
