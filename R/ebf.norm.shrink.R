# Calculate normal EBFs with shrinkage

ebf.norm.shrink <- function(x, s, index, xmin, xmax, points, pi0=0, complement=FALSE) {

  if (complement==FALSE)
    tails = is.infinite(xmin) + is.infinite(xmax)
  else
    tails = is.finite(xmin) + is.finite(xmax)

  # posterior marginal likelihood
  pml = rep(0, length(x))
  for(i in index) {

    # add current test onto subset of points
    points.save = points
    if (!(i %in% points)) points = c(i, points)

    # Laplace approximation
    ivw.var = 1/s[i]^2 + 1/s[points]^2
    ivw.mean = (x[i]/s[i]^2 + x[points]/s[points]^2) / ivw.var
    laplace = dnorm(x[i], ivw.mean, s[i]) *
      dnorm(x[points], ivw.mean, s[points]) *
      sqrt(2 * pi / ivw.var)

    area1 = pnorm(xmax, ivw.mean, 1/sqrt(ivw.var)) -
      pnorm(xmin, ivw.mean, 1/sqrt(ivw.var))
    area2 = pnorm(xmax, x[points], s[points]) - pnorm(xmin, x[points], s[points])
    if (complement == TRUE) {
      area1 = pnorm(xmin, ivw.mean, 1/sqrt(ivw.var)) +
        pnorm(xmax, ivw.mean, 1/sqrt(ivw.var), lower=F)
      area2 = pnorm(xmin, x[points], s[points]) + pnorm(xmax, x[points], s[points], lower=F)
    }

    ix = match(i,points)
    laplace[ix] = laplace[ix] / exp(0.5*tails/2)
    laplace[-ix] = laplace[-ix] * (1-pi0)
    area2[-ix] = area2[-ix] * (1-pi0)

    pml[i] = sum(laplace*area1) / sum(area2)

    points = points.save
  }
  pml[index]
}
