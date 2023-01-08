# Calculate binomial EBFs with shrinkage

ebf.binom.shrink <- function(x, size, index, xmin, xmax, shape, points, pi0=0,
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
        area = pbeta(xmax, x[i]+x[j]+shape, size[i]+size[j]-x[i]-x[j]+shape) -
          pbeta(xmin, x[i]+x[j]+shape, size[i]+size[j]-x[i]-x[j]+shape)
      else
        area = pbeta(xmin, x[i]+x[j]+shape, size[i]+size[j]-x[i]-x[j]+shape) +
          pbeta(xmax, x[i]+x[j]+shape, size[i]+size[j]-x[i]-x[j]+shape, lower=F)

      area1 = c(area1, exp(log(area) + lchoose(size[i], x[i]) +
                             lbeta(x[i]+x[j]+shape, size[i]+size[j]-x[i]-x[j]+shape) -
                             lbeta(x[j]+shape, size[j]-x[j]+shape)))
    }

    # normalising terms
    if (complement == FALSE)
      area2 = pbeta(xmax, x[points]+shape, size[points]-x[points]+shape) -
      pbeta(xmin, x[points]+shape, size[points]-x[points]+shape)
    else
      area2 = pbeta(xmin, x[points]+shape, size[points]-x[points]+shape) +
      pbeta(xmax, x[points]+shape, size[points]-x[points]+shape, lower=F)

    # adjust the diagonal elements for bias
    # bias
    if (xmax-xmin+complement==1 & shape==1) {
      if (size[i]<=length(binom.bias)) bias = binom.bias[size[i]]
      else bias = 0.5
    }
    else bias = compute.binom.bias(size[i], xmin, xmax, shape, complement)

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
