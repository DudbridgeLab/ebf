#' EBF bias for binomial tests
#'
#' Bias in the log posterior marginal likelihood for binomial tests,
#' used for calculating empirical Bayes factors.
#' Bias is computed by the \code{\link{compute.binom.bias}} function using enumeration.
#' Sample sizes range from 1 to 100.
#' The uniform prior is used.
#'
#' The bias depends on the true success probability.
#' As this is unknown, the bias is maximised over the probability to give
#' a conservative estimate.
#'
#' @format A vector with the bias for each sample size from 1 to 100.
#'
#' @template references
#'
"binom.bias"
