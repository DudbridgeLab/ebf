# Calculate simple Poisson EBFs

ebf.poisson.simple <- function(x, interval, xmin, xmax, complement=FALSE) {

  # posterior mean likelihood
  if (complement == FALSE) {
    area1 = (pgamma(xmax, 2*x, 2*interval) - pgamma(xmin, 2*x, 2*interval)) *
      gamma(2*x) / 4^x / gamma(x+1) / gamma(x)
    area2 = pgamma(xmax, x, interval) - pgamma(xmin, x, interval)
  } else {
    area1 = (pgamma(xmin, 2*x, 2*interval) + pgamma(xmax, 2*x, 2*interval, lower=F)) *
      gamma(2*x) / 4^x / gamma(x+1) / gamma(x)
    area2 = pgamma(xmin, x, interval) + pgamma(xmax, x, interval, lower=F)
  }

  # maximum bias corresponds to lambda=6.177394
  bias = 0.5406853

  # EBF
  area1 / area2 / exp(bias * area2)
}

