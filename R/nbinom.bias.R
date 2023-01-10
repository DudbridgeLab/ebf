#' EBF bias for negative binomial tests
#'
#' Bias in the log posterior marginal likelihood for negative binomial tests,
#' used for calculating empirical Bayes factors.
#' Bias is computed by the \code{\link{compute.nbinom.bias}} function using Monte Carlo simulation.
#' Minimum sample sizes (ie fixed number of events) range from 1 to 100.
#' The uniform prior is used.
#'
#' @format A vector with the bias for each sample size from 1 to 100.
#'
#' @template references
#'
"nbinom.bias"
