# Calculate t EBFs with shrinkage

ebf.t.shrink <- function(x, se, xmin, xmax, df, points, index, complement=FALSE) {

  if (FALSE) {
    # posterior mean likelihoods
    pml = matrix(nrow=length(x), ncol=length(x))
    for(i in 1:length(x))
      for(j in 1:length(x))
        if(complement == FALSE) {
          pml[i,j] = t.product.integral(x[i], se[i], x[j], se[j], df, xmin, xmax)
        } else {
          pml[i,j] = t.product.integral(x[i], se[i], x[j], se[j], df, -Inf, xmin) +
            t.product.integral(x[i], se[i], x[j], se[j], df, xmax, Inf)
        }

    # normalising terms
    if (complement == FALSE)
      area2 = mean(pt((xmax-x)/se, df) - pt((xmin-x)/se, df))
    else
      area2 = mean(pt((xmin-x)/se, df) + pt((xmax-x)/se, df, lower=F))

    # adjust the diagonal elements for bias
    if (df<=nrow(ebf::t.bias)) bias = ebf::t.bias$bias[df]
    else bias = 0.5
    for(i in 1:length(x))
      pml[i,i] = pml[i,i] / exp(area2 * bias)

    # form the EBFs
    apply(pml, 1, mean) / area2
  }

  if (TRUE) {
    pml = rep(0, length(index))
    for(i in 1:length(index)) {
      points.save = points
      if (!(index[i] %in% points)) points = c(index[i], points)

      t.integral = rep(0, length(points))

      for(j in 1:length(points)) {
        if(complement == FALSE) {
          t.integral[j] = t.product.integral(x[index[i]], se[index[i]],
                                             x[points[j]], se[points[j]], df, xmin, xmax)
        } else {
          t.integral[j] = t.product.integral(x[index[i]], se[index[i]],
                                             x[points[j]], se[points[j]], df, -Inf, xmin) +
            t.product.integral(x[index[i]], se[index[i]],
                               x[points[j]], se[points[j]], df, xmax, Inf)
        }
      }

      # normalising terms
      if (complement == FALSE)
        area2 = mean(pt((xmax-x[points])/se[points], df) -
                       pt((xmin-x[points])/se[points], df))
      else
        area2 = mean(pt((xmin-x[points])/se[points], df) +
                       pt((xmax-x[points])/se[points], df, lower=F))

      # adjust the diagonal elements for bias
      if (df<=nrow(ebf::t.bias)) bias = ebf::t.bias$bias[df]
      else bias = 0.5
      t.integral[match(index[i],points)] = t.integral[match(index[i],points)] / exp(area2 * bias)

      # form the EBFs
      pml[i] = mean(t.integral) / area2
      points = points.save
    }
    pml
  }

}
