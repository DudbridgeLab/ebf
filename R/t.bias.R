#' EBF bias for t tests
#'
#' Bias in the log posterior marginal likelihood for t tests,
#' used for calculating empirical Bayes factors.
#' Bias is computed by the \code{\link{compute.t.bias}} function using numerical integration.
#' Degrees of freedom range from 1 to 100.
#'
#' @format A vector with the bias for each degree of freedom from 1 to 100.
#'
#' @template references
"t.bias"
