# Calculate negative binomial EBFs with shrinkage

ebf.nbinom.shrink <- function(x, size, index, xmin, xmax, points, shape, complement=FALSE) {

  # posterior marginal likelihood
  pml = NULL
  for(i in index) {

    # add current test onto subset of points
    points.save = points
    if (!(i %in% points)) points = c(i, points)

        area1 = NULL
    for(j in points) {
      if (complement == FALSE)
        area = pbeta(xmax, size[i]+size[j]+shape, x[i]+x[j]+shape) -
          pbeta(xmin, size[i]+size[j]+shape, x[i]+x[j]+shape)
      else
        area = pbeta(xmin, size[i]+size[j]+shape, x[i]+x[j]+shape) +
          pbeta(xmax, size[i]+size[j]+shape, x[i]+x[j]+shape, lower=F)

      area1 = c(area1, exp(log(area) + lchoose(x[i]+size[i]-1, size[i]-1) +
                             lbeta(size[i]+size[j]+shape, x[i]+x[j]+shape) -
                             lbeta(size[j]+shape, x[j]+shape)))
    }

    # normalising terms
    if (complement == FALSE)
      area2 = mean(pbeta(xmax, size[points]+shape, x[points]+shape) -
                     pbeta(xmin, size[points]+shape, x[points]+shape))
    else
      area2 = mean(pbeta(xmin, size[points]+shape, x[points]+shape) +
                     pbeta(xmax, size[points]+shape, x[points]+shape, lower=F))

    # adjust the diagonal elements for bias
    if (size[i] <= length(nbinom.bias) & shape == 1)
      bias = ebf::nbinom.bias[size[i]]
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
