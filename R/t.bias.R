#' EBF bias for t tests
#'
#' Bias in the posterior mean likelihood for t tests,
#' used for calculating empirical Bayes factors.
#' Bias is computed by the \code{compute.t.bias} function using numerical integration.
#' Degrees of freedom range from 1 to 100.
#'
#' @format A data frame with degrees of freedom in the first column
#' and estimated bias in the second .
"t.bias"
