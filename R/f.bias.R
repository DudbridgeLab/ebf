#' EBF bias for F tests
#'
#' Bias in the log posterior marginal likelihood for F tests,
#' used for calculating empirical Bayes factors.
#' Bias is computed by the \code{\link{compute.f.bias}} function using numerical integration.
#' Degrees of freedom range from 1 to 100.
#'
#' @format A matrix with the bias for each degree of freedom from 1 to 100
#' in either numerator or denominator.
#'
#' @template references
"f.bias"
