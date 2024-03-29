% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/t.bias.R
\docType{data}
\name{t.bias}
\alias{t.bias}
\title{EBF bias for t tests}
\format{
A vector with the bias for each degree of freedom from 1 to 100.
}
\usage{
t.bias
}
\description{
Bias in the log posterior marginal likelihood for t tests,
used for calculating empirical Bayes factors.
Bias is computed by the \code{\link{compute.t.bias}} function using numerical integration.
Degrees of freedom range from 1 to 100.
}
\references{
Dudbridge F (2023) Empirical Bayes factors for common hypothesis tests.
arXiv:2301.11057
}
\author{
Frank Dudbridge
}
\keyword{datasets}
