# Calculate normal EBFs with shrinkage

ebf.norm.shrink <- function(x, se, xmin, xmax, index, complement=FALSE) {

  # IVW means of each pair of elements in x
  ivw.var = sapply(1/se^2,"+",1/se[index]^2)
  ivw.mean = sapply(x/se^2, "+", (x/se^2)[index]) / ivw.var

  # Laplace approximations for posterior mean likelihoods
  laplace = dnorm(x[index], ivw.mean, se[index]) *
    t(dnorm(x, t(ivw.mean), se)) * sqrt(2*pi/ivw.var)

  area1 = pnorm(xmax, ivw.mean, 1/sqrt(ivw.var)) -
    pnorm(xmin, ivw.mean, 1/sqrt(ivw.var))
  area2 = pnorm(xmax, x, se) - pnorm(xmin, x, se)
  if (complement == TRUE) {
    area1 = 1-area1
    area2 = 1-area2
    area1 = pnorm(xmin, ivw.mean, 1/sqrt(ivw.var)) +
      pnorm(xmax, ivw.mean, 1/sqrt(ivw.var), lower=F)
    area2 = pnorm(xmin, x, se) + pnorm(xmax, x, se, lower=F)

  }
  laplace = laplace * area1 / t(matrix(area2, nrow=length(x), ncol=nrow(laplace)))

  # adjust the diagonal elements for bias
  for(i in 1:length(index))
    laplace[i,index[i]] = laplace[i,index[i]] / exp(area2[i]/2)

    # form the EBFs
  apply(laplace, 1, mean)

}
