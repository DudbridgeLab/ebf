#' @param s Vector of standard deviations or standard errors associated with \code{x}.
#'
#' @param h0 If a scalar, the value of a point null hypothesis.
#' If a vector with two elements, the lower and upper
#' bounds of the null hypothesis.
#'
#' @param h1 If a scalar, the value of a point alternative hypothesis.
#' If a vector with two elements, the lower and upper bounds of
#' the alternative hypothesis.  If \code{NULL} (default), the alternative
#' hypothesis is the complement of \code{h0}.
#'
#' @details
#' The default test has \code{h0=0}, with two-sided alternative.
#' For one-sided alternatives, use \code{h1=c(0,Inf)} or \code{h1=c(-Inf,0)}.
#' To test positive values against negative values, use \code{h0=c(-Inf,0)},
#' in which case \code{h1} defaults to \code{c(0,Inf)}.
#' To test negative values against positive values, use \code{h0=c(0,Inf)}.

