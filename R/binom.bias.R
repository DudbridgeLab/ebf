#' EBF bias for binomial tests
#'
#' Bias in the posterior marginal likelihood for binomial tests,
#' used for calculating empirical Bayes factors.
#' Bias is computed by the \code{compute.binom.bias} function using enumeration.
#' Sample sizes range from 1 to 100.
#' The uniform prior is used.
#'
#' The bias depends on the true success probability.
#' As this is unknown, the maximum bias over success probabilities is taken
#' as a conservative estimate.
#'
#' @format A data frame with the bias for each sample size from 1 to 100.
#'
"binom.bias"
