# Calculate simple binomial EBFs

ebf.binom.simple <- function(x, size, xmin, xmax, shape, complement=FALSE) {

  # posterior mean likelihood
  area1 = rep(0, length(x))
  bias = rep(0, length(x))
  for(i in 1:length(x)) {
    if (complement == FALSE) {
      area1[i] = pbeta(xmax, 2*x[i]+shape, 2*(size[i]-x[i])+shape) -
        pbeta(xmin, 2*x[i]+shape, 2*(size[i]-x[i])+shape)
    } else {
      area1[i] = pbeta(xmin, 2*x[i]+shape, 2*(size[i]-x[i])+shape) +
        pbeta(xmax, 2*x[i]+shape, 2*(size[i]-x[i])+shape, lower=F)
    }
    area1[i] = exp(log(area1[i]) + lchoose(size[i], x[i]) +
                     lbeta(2*x[i]+shape, 2*(size[i]-x[i])+shape) -
                     lbeta(x[i]+shape, size[i]-x[i]+shape))
    # bias
    if (shape == 1 & size[i] <= length(binom.bias) &
        (complement == FALSE & xmin == 0 & xmax == 1 |
         complement == TRUE & xmin == xmax)) {
      bias[i] = binom.bias[size[i]]
    }
    else {
      bias[i] = compute.binom.bias(size[i], xmin, xmax, shape, complement)
    }
  }

  # normalising term
  if (complement == FALSE ) {
    area2 = pbeta(xmax, x+shape, (size-x)+shape) -
      pbeta(xmin, x+shape, (size-x)+shape)
  } else {
    area2 = pbeta(xmin, x+shape, (size-x)+shape) +
      pbeta(xmax, x+shape, (size-x)+shape, lower=F)
  }

  # EBF
  area1 / area2 / exp(bias)
}

