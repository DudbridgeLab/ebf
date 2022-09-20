# Calculate t EBFs with shrinkage

ebf.t.shrink <- function(x, se, xmin, xmax, df, index, complement=FALSE) {

  # posterior mean likelihoods
  pml = matrix(0, nrow=length(index), ncol=length(x))
  for(i in 1:length(index))
    for(j in 1:length(x))
      if(complement == FALSE) {
        pml[i,j] = t.product.integral(x[index[i]], se[index[i]],
                                      x[j], se[j], df, xmin, xmax)
      } else {
        pml[i,j] = t.product.integral(x[index[i]], se[index[i]],
                                      x[j], se[j], df, -Inf, xmin) +
          t.product.integral(x[index[i]], se[index[i]],
                             x[j], se[j], df, xmax, Inf)
      }

  # normalising terms
  if (complement == FALSE)
    area2 = pt((xmax-x)/se, df) - pt((xmin-x)/se, df)
  else
    area2 = pt((xmin-x)/se, df) + pt((xmax-x)/se, df, lower=F)
  pml = pml / t(matrix(area2, nrow=length(x), ncol=nrow(pml)))

  # adjust the diagonal elements for bias
  if (df<=nrow(ebf::t.bias)) bias = ebf::t.bias$bias[df]
  else bias = 0.5
  for(i in 1:length(index))
    pml[i,index[i]] = pml[i,index[i]] / exp(area2[index[i]] * bias)

  # form the EBFs
  apply(pml, 1, mean)

}
