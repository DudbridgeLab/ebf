# Calculate normal EBFs with shrinkage

ebf.norm.shrink <- function(x, se, xmin, xmax, points, index, complement=FALSE) {

  if (FALSE) {
    # IVW means of each pair of elements in x
    ivw.var = sapply(1/se[points]^2,"+",1/se^2)
    ivw.mean = sapply(x[points]/se[points]^2, "+", x/se^2) / ivw.var

    # Laplace approximations for posterior mean likelihoods
    laplace = dnorm(x, ivw.mean, se) *
      t(dnorm(x[points], t(ivw.mean), se[points])) * sqrt(2*pi/ivw.var)

    area1 = pnorm(xmax, ivw.mean, 1/sqrt(ivw.var)) -
      pnorm(xmin, ivw.mean, 1/sqrt(ivw.var))
    area2 = mean(pnorm(xmax, x[points], se[points]) - pnorm(xmin, x[points], se[points]))
    if (complement == TRUE) {
      area1 = pnorm(xmin, ivw.mean, 1/sqrt(ivw.var)) +
        pnorm(xmax, ivw.mean, 1/sqrt(ivw.var), lower=F)
      area2 = mean(pnorm(xmin, x[points], se[points]) + pnorm(xmax, x[points], se[points], lower=F))

    }
    laplace = laplace * area1

    # adjust the diagonal elements for bias
    for(i in 1:length(x))
      if (i %in% points)
        laplace[i, match(i,points)] = laplace[i, match(i,points)] / exp(area2/2)

    # form the EBFs
    apply(laplace, 1, mean) / area2
  }

  if (TRUE) {

    # posterior mean likelihood
    #print(index)
    pml = rep(0, length(index))
    for(i in 1:length(index)) {
      points.save = points
      if (!(index[i] %in% points)) points = c(index[i], points)

      ivw.var = 1/se[index[i]]^2 + 1/se[points]^2
      ivw.mean = (x[index[i]]/se[index[i]]^2 + x[points]/se[points]^2) / ivw.var

      laplace = dnorm(x[index[i]], ivw.mean, se[index[i]]) *
        dnorm(x[points], ivw.mean, se[points]) *
        sqrt(2 * pi / ivw.var)

      area1 = pnorm(xmax, ivw.mean, 1/sqrt(ivw.var)) -
        pnorm(xmin, ivw.mean, 1/sqrt(ivw.var))
      area2 = mean(pnorm(xmax, x[points], se[points]) - pnorm(xmin, x[points], se[points]))
      if (complement == TRUE) {
        area1 = pnorm(xmin, ivw.mean, 1/sqrt(ivw.var)) +
          pnorm(xmax, ivw.mean, 1/sqrt(ivw.var), lower=F)
        area2 = mean(pnorm(xmin, x[points], se[points]) + pnorm(xmax, x[points], se[points], lower=F))
      }

      laplace[match(index[i],points)] = laplace[match(index[i],points)] / exp(area2*0.5)

      pml[i] = mean(laplace * area1) / area2

      points = points.save
    }
    pml
  }
}
