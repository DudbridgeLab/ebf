# Calculate Poisson EBFs with shrinkage

ebf.poisson.shrink <- function(x, interval, index, xmin, xmax, points, complement=FALSE) {

  # posterior marginal likelihood
  pml = NULL
  for(i in index) {

    # add current test onto subset of points
    points.save = points
    if (!(i %in% points)) points = c(i, points)

    area1 = NULL
    for(j in points) {
      if (complement == FALSE)
        area = integrate(function(lambda) {
          dpois(x[i], lambda*interval[i]) * dgamma(lambda, x[j]/interval[j], 1)
        }, xmin, xmax)$value
      else
        area = integrate(function(lambda) {
          dpois(x[i], lambda*interval[i]) * dgamma(lambda, x[j]/interval[j], 1)
        }, 0, xmin)$value +
          integrate(function(lambda) {
            dpois(x[i], lambda*interval[i]) * dgamma(lambda, x[j]/interval[j], 1)
          }, xmax, Inf)$value

      area1 = c(area1, area)
    }

    # normalising terms
    if (complement == FALSE) {
      area2 = mean(pgamma(xmax, x[points]/interval[points], 1) -
                      pgamma(xmin, x[points]/interval[points], 1) )
    } else {
      area2 = mean(pgamma(xmin, x[i]/interval[i], 1) +
                     pgamma(xmax, x[i]/interval[i], 1, lower=F))
    }

    # adjust the diagonal elements for bias
    # maximum bias corresponds to lambda=7.301
    bias = 0.5346

    area1[match(i,points)] = area1[match(i,points)] /
      exp(area2 * bias)

    # form the EBFs
    pml = c(pml, mean(area1) / area2)

    # restore subset of points
    points = points.save
  }

  pml
}
