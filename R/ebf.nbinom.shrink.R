# Calculate negative binomial EBFs with shrinkage

ebf.nbinom.shrink <- function(x, size, index, xmin, xmax, shape, points, pi0=0,
                              complement=FALSE) {

  # posterior marginal likelihood
  pml = rep(0, length(x))
  for(i in index) {

    # add current test onto subset of points
    points.save = points
    if (!(i %in% points)) points = c(i, points)

    area1 = NULL
    for(j in points) {
      if (complement == FALSE)
        area = pbeta(xmax, size[i]+size[j]+0*shape, x[i]+x[j]+shape) -
          pbeta(xmin, size[i]+size[j]+0*shape, x[i]+x[j]+shape)
      else
        area = pbeta(xmin, size[i]+size[j]+0*shape, x[i]+x[j]+shape) +
          pbeta(xmax, size[i]+size[j]+0*shape, x[i]+x[j]+shape, lower=F)

      area1 = c(area1, exp(log(area) + lchoose(x[i]+size[i]-1, size[i]-1) +
                             lbeta(size[i]+size[j]+0*shape, x[i]+x[j]+shape) -
                             lbeta(size[j]+0*shape, x[j]+shape)))
    }

    # normalising terms
    if (complement == FALSE)
      area2 = pbeta(xmax, size[points]+0*shape, x[points]+shape) -
      pbeta(xmin, size[points]+0*shape, x[points]+shape)
    else
      area2 = pbeta(xmin, size[points]+0*shape, x[points]+shape) +
      pbeta(xmax, size[points]+0*shape, x[points]+shape, lower=F)

    # adjust the diagonal elements for bias
    if (size[i] <= length(nbinom.bias) & shape == 1)
      bias = ebf::nbinom.bias[size[i]]
    else
      bias = 0.5

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
