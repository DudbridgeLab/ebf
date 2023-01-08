# Calculate simple Poisson EBFs

ebf.poisson.simple <- function(x, interval, xmin, xmax, shape, rate, complement=FALSE) {

  # posterior mean likelihood
  if (complement == FALSE) {
    area1 = (pgamma(xmax, 2*x+shape, 2*interval+rate) -
               pgamma(xmin, 2*x+shape, 2*interval+rate)) *
      exp(lgamma(2*x+shape) + (x+shape)*log(interval+rate) -
            lfactorial(x) - lgamma(x+shape) - (2*x+shape)*log(2*interval+rate))
    area2 = pgamma(xmax, x+shape, interval+rate) - pgamma(xmin, x+shape, interval+rate)
  } else {
    area1 = (pgamma(xmin, 2*x+shape, 2*interval+rate) +
               pgamma(xmax, 2*x+shape, 2*interval+rate, lower=F)) *
      exp(lgamma(2*x+shape) + (x+shape)*log(interval+rate) -
            lfactorial(x) - lgamma(x+shape) - (2*x+shape)*log(2*interval+rate))
    area2 = pgamma(xmin, x+shape, interval+rate) +
      pgamma(xmax, x+shape, interval+rate, lower=F)
  }

  if (is.finite(xmax) == complement) bias = 0.5
  else bias = 0

  # EBF
  area1 / area2 / exp(bias)
}

