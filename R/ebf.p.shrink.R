# Calculate EBFs for p-values with shrinkage to non-parametric posterior

ebf.p.shrink <- function(x, index, points, pi0=0) {

  # transform to gamma parameter
  q = -log(1-x)

  # posterior marginal likelihood
  pml = rep(0, length(x))

  for(i in index) {

    # add current test onto subset of points
    points.save = points
    if (!(i %in% points)) points = c(i, points)

    # integral of beta(1,b)*gamma(b,2,-log(1-p))
    area1 = exp(2*log(q[points]) + log(2) + q[i] - 3*log(q[i]+q[points])) *
      pgamma(1, 3, q[i]+q[points], lower=F)
    area2 = pgamma(1, 2, q[points], lower=F)

    ix = match(i,points)
    area1[ix] = area1[ix] / (5/2) / (1-pi0)
    area2[ix] = area2[ix] / (1-pi0)

    # form the EBFs
    pml[i] = sum(area1) / sum(area2)

    # restore subset of points
    points = points.save
  }

  pml
}
