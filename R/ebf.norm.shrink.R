# Calculate normal EBFs with shrinkage

ebf.norm.shrink <- function(x, se, index, xmin, xmax, points, complement=FALSE) {

  # posterior marginal likelihood
  pml = NULL
  for(i in index) {

    # add current test onto subset of points
    points.save = points
    if (!(i %in% points)) points = c(i, points)

    # Laplace approximation
    ivw.var = 1/se[i]^2 + 1/se[points]^2
    ivw.mean = (x[i]/se[i]^2 + x[points]/se[points]^2) / ivw.var

    laplace = dnorm(x[i], ivw.mean, se[i]) *
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

    laplace[match(i,points)] = laplace[match(i,points)] / exp(area2*0.5)

    pml = c(pml, mean(laplace * area1) / area2)

    points = points.save
  }
  pml
}

