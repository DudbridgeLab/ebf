# Calculate simple Poisson EBFs

ebf.poisson.simple <- function(x, interval, xmin, xmax, complement=FALSE) {

  # posterior mean likelihood
  area1 = rep(0, length(x))
  bias = rep(0, length(x))
  for(i in 1:length(x)) {
    if (complement == FALSE) {
      area1[i] = integrate(function(lambda) {
        dpois(x[i], lambda*interval[i]) * dgamma(lambda, x[i]/interval[i], 1)
      }, xmin, xmax)$value
    } else {
      area1[i] = integrate(function(lambda) {
        dpois(x[i], lambda*interval[i]) * dgamma(lambda, x[i]/interval[i], 1)
      }, 0,  xmin)$value +
        integrate(function(lambda) {
          dpois(x[i], lambda*interval[i]) * dgamma(lambda, x[i]/interval[i], 1)
        }, xmax, Inf)$value
    }

  }

  # maximum bias corresponds to lambda=7.301
  bias = 0.5346

  # normalising term
  if (complement == FALSE ) {
    area2 = pgamma(xmax, x[i]/interval[i], 1) -
      pgamma(xmin, x[i]/interval[i], 1)
  } else {
    area2 = pgamma(xmin, x[i]/interval[i], 1) +
      pgamma(xmax, x[i]/interval[i], 1, lower=F)
  }

  # EBF
  area1 / area2 / exp(bias * area2)
}

