# Calculate simple binomial EBFs

ebf.binom.simple <- function(x, size, xmin, xmax, shape, complement=FALSE) {

  # posterior mean likelihood
  if (complement == FALSE) {
    area1 = pbeta(xmax, 2*x+shape, 2*(size-x)+shape) -
      pbeta(xmin, 2*x+shape, 2*(size-x)+shape)
  } else {
    area1 = pbeta(xmin, 2*x+shape, 2*(size-x)+shape) +
      pbeta(xmax, 2*x+shape, 2*(size-x)+shape, lower=F)
  }
  area1 = exp(log(area1) + lchoose(size, x) +
                lbeta(2*x+shape, 2*(size-x)+shape) -
                lbeta(x+shape, size-x+shape))

  # normalising term
  if (complement == FALSE ) {
    area2 = pbeta(xmax, x+shape, (size-x)+shape) -
      pbeta(xmin, x+shape, (size-x)+shape)
  } else {
    area2 = pbeta(xmin, x+shape, (size-x)+shape) +
      pbeta(xmax, x+shape, (size-x)+shape, lower=F)
  }

  # bias
  bias = rep(0,length(x))
  for(i in 1:length(x)) {
    if (xmax-xmin+complement==1 & shape==1) {
      if (size[i]<=length(binom.bias)) bias[i] = binom.bias[size[i]]
      else bias[i] = 0.5
    }
    else bias[i] = compute.binom.bias(size[i], xmin, xmax, shape, complement)
  }

  # EBF
  area1 / area2 / exp(bias)
}

