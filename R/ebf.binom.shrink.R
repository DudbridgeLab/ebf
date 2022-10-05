# Calculate binomial EBFs with shrinkage

ebf.binom.shrink <- function(x, size, index, xmin, xmax, points, shape, complement=FALSE) {

  # posterior marginal likelihood
  pml = NULL
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
      area2 = mean(pbeta(xmax, x[points]+shape, size[points]-x[points]+shape) -
                     pbeta(xmin, x[points]+shape, size[points]-x[points]+shape))
    else
      area2 = mean(pbeta(xmin, x[points]+shape, size[points]-x[points]+shape) +
                     pbeta(xmax, x[points]+shape, size[points]-x[points]+shape, lower=F))

    # adjust the diagonal elements for bias
    if (size[i] <= length(binom.bias) & shape == 1)
      bias = ebf::binom.bias[size[i]]
    else
      bias = 0.5

    area1[match(i,points)] = area1[match(i,points)] /
      exp(area2 * bias)

    # form the EBFs
    pml = c(pml, mean(area1) / area2)

    # restore subset of points
    points = points.save
  }

  pml
}
