# Calculate t EBFs with shrinkage

ebf.t.shrink <- function(x, s, df, index, xmin, xmax, points, pi0=0, complement=FALSE) {

  # posterior marginal likelihood
  pml = rep(0, length(x))
  for(i in index) {

    # add current test onto subset of points
    points.save = points
    if (!(i %in% points)) points = c(i, points)

    area1 = NULL
    for(j in points) {
      if(complement == FALSE) {
        area1 = c(area1, t.product.integral(x[i], s[i], df[i],
                                            x[j], s[j], df[j], xmin, xmax))
      } else {
        area1 = c(area1, t.product.integral(x[i], s[i], df[i],
                                            x[j], s[j], df[j], -Inf, xmin) +
                    t.product.integral(x[i], s[i], df[i],
                                       x[j], s[j], df[j], xmax, Inf))
      }
    }

    # normalising terms
    if (complement == FALSE)
      area2 = pt((xmax-x[points])/s[points], df[points]) -
      pt((xmin-x[points])/s[points], df[points])
    else
      area2 = pt((xmin-x[points])/s[points], df[points]) +
      pt((xmax-x[points])/s[points], df[points], lower=F)

    # adjust the diagonal elements for bias
    if (df[i] <= length(t.bias))
      bias = ebf::t.bias[df[i]]
    else
      bias = 0.5

    ix = match(i,points)
    area1[ix] = area1[ix] / exp(bias)
    area1[-ix] = area1[-ix] * (1-pi0)
    area2[-ix] = area2[-ix] * (1-pi0)

    # form the EBFs
    pml[i] = sum(area1) / sum(area2)

    # restore subset of points
    points = points.save
  }

  pml
}
