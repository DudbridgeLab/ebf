# Calculate Poisson EBFs with shrinkage

ebf.poisson.shrink <- function(x, interval, index, xmin, xmax, shape, rate,
                               points, pi0=0, complement=FALSE) {

  # posterior marginal likelihood
  pml = rep(0, length(x))

  #bias
  if (is.finite(xmax) == complement) bias = 0.5
  else bias = 0

  for(i in index) {

    # add current test onto subset of points
    points.save = points
    if (!(i %in% points)) points = c(i, points)

    if (complement == FALSE) {
      area1 = dnbinom(x[i], x[points]+shape, interval[i]/(interval[i]+interval[points]+rate)) *
        (pgamma(xmax, x[i]+x[points]+shape, interval[i]+interval[points]+rate) -
           pgamma(xmin, x[i]+x[points]+shape, interval[i]+interval[points]+rate))
      area2 = pgamma(xmax, x[points]+shape, interval[points]+rate) -
        pgamma(xmin, x[points]+shape, interval[points]+rate)
    }
    else {
      area1 = dnbinom(x[i], x[points]+shape, interval[i]/(interval[i]+interval[points]+rate)) *
        (pgamma(xmin, x[i]+x[points]+shape, interval[i]+interval[points]+rate) +
           pgamma(xmax, x[i]+x[points]+shape, interval[i]+interval[points]+rate, lower=F))
      area2 = pgamma(xmin, x[points]+shape, interval[points]+rate) +
        pgamma(xmax, x[points]+shape, interval[points]+rate, lower=F)
    }

    # adjust the diagonal elements for bias
    ix = match(i,points)
    area1[ix] = area1[ix] / exp(bias) / (1-pi0)
    area2[ix] = area2[ix] / (1-pi0)

    # form the EBFs
    pml[i] = sum(area1) / sum(area2)

    # restore subset of points
    points = points.save
  }

  pml
}
