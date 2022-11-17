# Calculate simple Poisson EBFs

ebf.poisson.simple <- function(x, interval, xmin, xmax, complement=FALSE) {

  # posterior mean likelihood
  if (complement == FALSE) {
    area1 = (pgamma(xmax, 2*x+1, 2*interval) - pgamma(xmin, 2*x+1, 2*interval)) *
      exp(lfactorial(2*x) - (2*x+1)*log(2) - lfactorial(x) - lfactorial(x))
    area2 = pgamma(xmax, x+1, interval) - pgamma(xmin, x+1, interval)
  } else {
    area1 = (pgamma(xmin, 2*x+1, 2*interval) + pgamma(xmax, 2*x+1, 2*interval, lower=F)) *
      exp(lfactorial(2*x) - (2*x+1)*log(2) - lfactorial(x) - lfactorial(x))
    area2 = pgamma(xmin, x+1, interval) + pgamma(xmax, x+1, interval, lower=F)
  }

  # maximum bias corresponds to lambda=6.411577
  bias = 0.5003284

  # EBF
  area1 / area2 / exp(bias * area2)
}

