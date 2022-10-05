# Calculate t EBFs with shrinkage

ebf.t.shrink <- function(x, se, df, index, xmin, xmax, points, complement=FALSE) {

  # posterior marginal likelihood
  pml = NULL
  for(i in index) {

    # add current test onto subset of points
    points.save = points
    if (!(i %in% points)) points = c(i, points)

    area1 = NULL
    for(j in points) {
      if(complement == FALSE) {
        area1 = c(area1, t.product.integral(x[i], se[i], df[i],
                                            x[j], se[j], df[j], xmin, xmax))
      } else {
        area1 = c(area1, t.product.integral(x[i], se[i], df[i],
                                            x[j], se[j], df[j], -Inf, xmin) +
                    t.product.integral(x[i], se[i], df[i],
                                       x[j], se[j], df[j], xmax, Inf))
      }
    }

    # normalising terms
    if (complement == FALSE)
      area2 = mean(pt((xmax-x[points])/se[points], df[points]) -
                     pt((xmin-x[points])/se[points], df[points]))
    else
      area2 = mean(pt((xmin-x[points])/se[points], df[points]) +
                     pt((xmax-x[points])/se[points], df[points], lower=F))

    # adjust the diagonal elements for bias
    if (df[i] <= length(t.bias))
      bias = ebf::t.bias[df[i]]
    else
      bias = 0.5

    area1[match(i,points)] = area1[match(i,points)] / exp(area2 * bias)

    # form the EBFs
    pml = c(pml, mean(area1) / area2)

    # restore subset of points
    points = points.save
  }

  pml
}
